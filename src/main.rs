use std::io;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use std::error::Error;
use std::env;

#[derive(Debug)]
struct FastqRecord {
    header: String,
    seq: String,
    qual: String,
}

impl FastqRecord {
    fn init() -> FastqRecord {
        FastqRecord {
            header: String::new(),
            seq: String::new(),
            qual: String::new(),
        }
    }
}

fn read_fastq_one(path_str: &str) -> Vec<FastqRecord> {
    let path = Path::new(path_str);
    //let path = Path::new(r"C:\Users\naupio\Desktop\mywork\testread.txt");
    let display = path.display();
    //let new_path = path.join("L858R.txt");

    //match path.to_str() {
    //    Some(s) => println!("{:?}", s),
    //    None => {},
    //}

    let mut file = match File::open(&path) {
        // The `description` method of `io::Error` returns a string that
        // describes the error
        Err(why) => panic!("couldn't open {}: {}", display,
                                                   Error::description(&why)),
        Ok(file) => file,
    };

    let mut s = String::new();
    match file.read_to_string(&mut s) {
        Err(why) => panic!("couldn't read {}: {}", display,
                                                   Error::description(&why)),
        Ok(_) => {},//print!("{} contains:\n{}", display, s),
    }

    let mut iter = s.lines();
    let mut cnt = 0;
    let mut fq_rd = FastqRecord::init();
    let mut fq_rd_vec:Vec<FastqRecord> = vec![];

    loop{
        //if cnt <= 4 {cnt += 1;}
        //else { cnt = 0; fq_rd = FastqRecord::init(); break;}
        match iter.next() {
            Some(it) => {
                match ( cnt % 4) {
                    0 => fq_rd.header = it.to_string(),
                    1 => fq_rd.seq = it.to_string(),
                    3 => fq_rd.qual = it.to_string(),

                    2 => {},
                    _ => {},
                } 
            }
            None => {break;},
        }

        if cnt % 4 == 3 {
            //println!("{:?}", fq_rd);
            fq_rd_vec.push(fq_rd);
            fq_rd = FastqRecord::init();
            //break;
        }

        //if cnt > 500 {println!("the line in read_fastq_one function havn't comment!"); break;}
        cnt += 1;
    }
    fq_rd_vec
}

fn read_fastq_pair(path_str1: &str, path_str2: &str) -> (Vec<FastqRecord>,Vec<FastqRecord>) {
    (read_fastq_one(path_str1),read_fastq_one(path_str2))
}

#[derive(Debug)]
struct SnpRecord {
    label: String,
    seq: String,
}

impl SnpRecord {
    fn init() -> SnpRecord {
        SnpRecord {
            label: String::new(),
            seq: String::new(),
        }
    }
}

fn read_snp_file(path_str: &str) -> Vec<SnpRecord> {
    // Create a path to the desired file
    let path = Path::new(path_str);
    let display = path.display();

    // Open the path in read-only mode, returns `io::Result<File>`
    let mut file = match File::open(&path) {
        // The `description` method of `io::Error` returns a string that
        // describes the error
        Err(why) => panic!("couldn't open {}: {}", display,
                                                   Error::description(&why)),
        Ok(file) => file,
    };
    // Read the file contents into a string, returns `io::Result<usize>`
    let mut s = String::new();
    match file.read_to_string(&mut s) {
        Err(why) => panic!("couldn't read {}: {}", display,
                                                   Error::description(&why)),
        Ok(_) => {}, //print!("{} contains:\n{}", display, s),
    }

    let mut iters = s.lines();
    let mut snp_rd = SnpRecord::init();
    let mut snp_rd_vec:Vec<SnpRecord> = vec![];
    let mut cnt = 0;

    loop {
        match iters.next() {
            Some(it) => {
               // println!("{:?}\n", it);
                match (cnt % 2) {
                    0 => snp_rd.label = it.to_string(),
                    1 => snp_rd.seq = it.to_string(),
                    _ => {},
                }

            },
            None => break,
        }
        if cnt % 2 == 1 {
            //println!("{:?}\n",snp_rd );
            snp_rd_vec.push(snp_rd);
            snp_rd = SnpRecord::init();
        }
        cnt += 1;
    }
    snp_rd_vec
}

fn edit_distance(ref_str: &str, tar_str: &str) -> usize {

    let len_ref = ref_str.chars().count();
    let len_tar = tar_str.chars().count();

    let row: Vec<usize> = vec![0; len_tar + 1];
    let mut matrix: Vec<Vec<usize>> = vec![row; len_ref + 1];

    // initialize string len_ref
    for i in 0..len_ref {
        matrix[i+1][0] = matrix[i][0] + 1;
    }

    // initialize string len_tar
    for i in 0..len_tar {
        matrix[0][i+1] = matrix[0][i] + 1;
    }

    // calculate matrix
    for (i, refc) in ref_str.chars().enumerate() {
        for (j, tarc) in tar_str.chars().enumerate() {
            let alternatives = [
                // deletion
                matrix[i][j+1] + 1,
                // insertion
                matrix[i+1][j] + 1,
                // match or substitution
                matrix[i][j] + if refc == tarc { 0 } else { 1 }];
            matrix[i+1][j+1] = *alternatives.iter().min().unwrap();
        }
    }

    matrix[len_ref][len_tar]
}

fn hamming_distance(s1: &str, s2: &str) -> u32{
    //println!("start hamming_distance(s1, s2)");
    if s1.len() != s2.len() {
        panic!("Both the length of the string is not equal.");
    }
    else {
        //println!("start else arm");
        let mut dis : u32 = 0;
        for i in 0..s1.len(){
            if s1.chars().nth(i) != s2.chars().nth(i) {
                dis += 1;
                //println!("{:?}",dis );
            }
        }
        dis
    }
}

fn fusion_string(r1_str: &str, r2_str: &str, mut fusion_num: usize) -> (bool, String) {
    let (flag, score) = overlap(r1_str, r2_str);
    if fusion_num < 4 {fusion_num = 4;}
    
    let r2_str = rev_str(r2_str);
    let r1_tail: &str;
    let r2_head: &str;

    let mut fu_s = String::new();
    if score >= fusion_num {
        r1_tail = &r1_str[r1_str.len()-flag..];
        r2_head = &r2_str[..flag];
        fu_s = process_n(r1_tail,r2_head);
    }
    else {
        return (false, fu_s)
    }

    let r1_head = &r1_str[..r1_str.len()-flag];
    let r2_tail = &r2_str[flag..];

    fu_s = r1_head.clone().to_string() + &fu_s + r2_tail;
    (true, fu_s)
}

fn process_n(s1: &str, s2: &str) ->String {
    if s1.len() != s2.len() {panic!("Both length of string is not equal!");}
    let mut ps = String::new();
    let mut s1 = s1.chars();
    let mut s2 = s2.chars();
    loop {
        let s2_ch = s2.next();
        let i = s1.next();
        match i{
            Some('N') => {
                match s2_ch {
                    Some(j) => {ps.push(j)},
                    None => break,
                }
            },
            Some(i) => {ps.push(i)},
            None => break,
        }
    }
    ps
}

fn overlap(r1_str: &str, r2_str: &str) ->(usize, usize) {
    let r1_string = r1_str.to_string(); 
    let r2_string = rev_str(r2_str);
    //println!("{:?}{:?}", r1_string, r2_string);

    match r1_str.len()==r2_str.len(){
        true => {},
        false => {return (0, 0)},
    }
    //println!("{:?}", pop_string(r1_string,r1_str.len()-1));

    //println!("start:");
    
    let mut flag = 0;
    let mut score = 0;
    
    for num in 0..r1_str.len()+1 {
        let ps_r1 = pop_string(r1_string.clone(), num);
        let ps_r1 = rev_str(&ps_r1);
        let ps_r2 = &r2_string[0..ps_r1.len()];
        //println!("{:?}{:?}", ps_r1, ps_r2);

        let (yes_no, f, sc) = is_similar(&ps_r1, &ps_r2);
        match yes_no {
            false => continue,
            true => {
                if sc >= score{
                    score = sc;
                    flag = f; 
                }
                else{ continue; } 
            },
        }

    }
    //println!("end!!");
    //println!("{:?}", (flag, score) );
    (flag, score)
}

fn pop_string(mut s: String, num: usize) -> String {
    let mut ps = String::new();
    for _ in 0..num {
        match s.pop() {
            Some(ch) => ps.push(ch),
            None =>  {println!("string is empty")},
        }
    }
    ps
}


fn is_similar(s1: &str, s2: &str)-> (bool, usize, usize) {
    let mut s1_chars = s1.chars();
    let mut s2_chars = s2.chars();
    {
        if s1.len() != s2.len() {
            panic!("Both the length of the string is not equal.");
        }
    }

    let mut flag = s1.len();
    let mut score = 0;
    for _ in 0..(s1.len()) {
        let next_s1ch = s1_chars.next();
        let next_s2ch = s2_chars.next();
        if next_s1ch != next_s2ch{ 
                if next_s1ch == Some('N') || next_s2ch == Some('N') {
                    continue;
                } 
                else{
                    flag = 0;break;    
                }
                
            }
        else {
            match next_s1ch {
                Some('N') => continue,
                _ => score += 1,
            }
        }
    }

    match flag{
        0 => (false, flag, 0),
        _ => (true, flag, score),
    }
}

fn rev_str( s: &str) ->String {
    let mut s1 = String::new();
    let mut s0 = s.to_string();
    let mut i = 0;
    while i < s.len() {
        match s0.pop(){
            Some(ch) => s1.push(ch),
            None => {},
        }
        i+=1;
    }
    //println!("{:?}", s1);
    s1
}

fn complement_str(s: &str) -> String {
    let mut s = s.chars();
    let mut cs = String::new();
    loop {
        let ch = match s.next() {
            Some('A') => 'T',
            Some('G') => 'C',
            Some('C') => 'G',
            Some('T') => 'A',
            Some(_) => 'N',
            None => break,
        };
        cs.push(ch);
    }
    cs
}

#[derive(Debug)]
struct SnpResult {
    snp_label: String,
    snp_seq: String,
    fastq_header: String,
    fastq_seq: String,
    fastq_qual: String,
}

impl SnpResult {
    fn init() -> SnpResult {
        SnpResult {
            snp_label: String::new(),
            snp_seq: String::new(),
            fastq_header: String::new(),
            fastq_seq: String::new(),
            fastq_qual: String::new(),
        }
    }
}

fn search_snp(it: &FastqRecord, snp: &SnpRecord, threshold_value: usize) -> (bool, SnpResult) {
    let len1 = it.seq.len();
    let len2 = snp.seq.len();
    //if(len1 <= len2) {continue;}
    //println!("len1 = {:?}, len2 = {:?}", len1, len2);
    let mut ed_dis: usize = len1+len2;
    for i in 0..len1-len2+1 {
        let temp_dis = edit_distance(&it.seq[i..i+len2], &snp.seq);
        //println!("{:?}\n", temp_dis);
        if temp_dis <= ed_dis {ed_dis = temp_dis;}
        else {continue;}
    }

    if(ed_dis <= threshold_value) {
        println!("edit_distance => {:?}\n", ed_dis);
        let snp_rs = SnpResult {
            snp_label: snp.label.clone(),
            snp_seq: snp.seq.clone(),
            fastq_header: it.header.clone(),
            fastq_seq: it.seq.clone(),
            fastq_qual: it.qual.clone(),
        };
        println!("{:?}\n", snp_rs);
        (true, snp_rs)
    }
    else {
        let snp_rs = SnpResult::init(); 
        (false, snp_rs) 
    }
}

fn run_with_one(path_str: &str, snp_path: &str, threshold_value: usize) -> Vec<SnpResult> {
    let r1_vec = read_fastq_one(path_str);
    let mut result_vec: Vec<SnpResult> = vec![];
    let snp_vec = read_snp_file(snp_path);

    for snp in snp_vec {
        let mut r1_vec_iter = r1_vec.iter();
        loop {
            match r1_vec_iter.next() {
                Some(it) => {
                    let len1 = it.seq.len();
                    let len2 = snp.seq.len();
                    if(len1 <= len2) {continue;}
                    let (flag, snp_rs) = search_snp(&it, &snp, threshold_value);
                    if flag{
                        result_vec.push(snp_rs);
                    }
                    else {
                        continue;
                    }
                    
                }
                None => break,
            }
        }
    }
    result_vec
}

fn run_with_pair(path_str1: &str, path_str2: &str, snp_path: &str, threshold_value: usize, fusion_num: usize) -> Vec<SnpResult> {
    let (r1_vec,r2_vec) = read_fastq_pair(path_str1, path_str2);
    let mut result_vec: Vec<SnpResult> = vec![];
    let snp_vec = read_snp_file(snp_path);

    for snp in snp_vec {
        let mut r1_vec_iter = r1_vec.iter();
        let mut r2_vec_iter = r2_vec.iter();
        loop {
            match (r1_vec_iter.next(), r2_vec_iter.next() ) {
                (Some(it1),Some(it2)) => {
                    //println!("test:\n");
                    let (isfu, fu_s) = fusion_string(&it1.seq, &it2.seq, fusion_num);
                    match isfu {
                        true => {
                            //println!("isfu:");
                            //println!("{:?}", fu_s);
                            let mut it = FastqRecord::init();
                            it.header = it1.header.clone() + &it2.header;
                            it.seq = fu_s;
                            //println!("{:?}", it);
                            let len1 = it.seq.len();
                            let len2 = snp.seq.len();
                            if(len1 <= len2) {continue;}

                            let (flag, snp_rs) = search_snp(&it, &snp, threshold_value);
                            if flag{
                                result_vec.push(snp_rs);
                            }
                            else {
                                continue;
                            }
                        }
                        false => {
                            let len1 = it1.seq.len();
                            let len2 = snp.seq.len();
                            let len3 = it2.seq.len();

                            if(len1 <= len2){continue;}
                            if(len3 <= len2){continue;}
                            let (flag1, snp_rs1) = search_snp(&it1, &snp, threshold_value);
                            let (flag2, snp_rs2) = search_snp(&it2, &snp, threshold_value);
                            //println!("{:?}======{:?}",flag1,flag2 );
                            match flag1 {
                                true => result_vec.push(snp_rs1),
                                false => {},
                            }
                            match flag2 {
                                true => result_vec.push(snp_rs2),
                                false => {continue;},
                            }
                        }
                    }
                }
                (Some(_), None) => {panic!("R2 fastq file is not long enough.");},
                (None, Some(_) ) => {panic!("R1 fastq file is not long enough.");},
                _ => {break;}
            }
        }
    }
    result_vec
}

fn main() {
    let args: Vec<String> = env::args().collect();    
    println!("{:?}", args);
    match args.len() >= 4 {
        true => {},
        false => panic!("the length of arguments vec is not bigger than 4"),
    }
    let ref r1 = args[1];
    let ref r2 = args[2];
    let ref snp = args[3];

    let path_str1 = r1; //  ./R1.fastq";
    let path_str2 = r2; //  ./R2.fastq";
    let snp_path = snp; //  ./L858R.txt";
    let result_vec = run_with_pair(&path_str1, &path_str2, &snp_path, 2, 4);
    for rv in result_vec {
        println!("{:?}", rv);
    }
}
