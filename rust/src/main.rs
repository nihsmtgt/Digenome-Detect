use std::process;
use rust_htslib::bam;
use rust_htslib::bam::Read;
// use rust_htslib::prelude::*;
// for cleavage score
use std::collections::VecDeque;
use std::io::Write;
use std::io::stdout;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
// detects digenome-seq breakpoint

fn main() {
    let margin: u32 = 150;
    let read_length: u32  = 150;
    let capacity: usize = read_length as usize * 2 + margin as usize * 2;
    let argv: Vec<String> = std::env::args().collect();
    let mut tname = String::from(&argv[1]);
    let mqfilter: u8 = String::from(&argv[3]).parse().unwrap();

    let mut bam = bam::IndexedReader::from_path(&argv[2]).unwrap();
    let header = bam.header();

    // for cleavage score
    let mut vec_depth: VecDeque<u32> = VecDeque::with_capacity(capacity);
    let mut vec_pos: VecDeque<u32> = VecDeque::with_capacity(capacity);
    let mut vec_forward_heads: VecDeque<u32> = VecDeque::with_capacity(capacity);
    let mut vec_reverse_heads: VecDeque<u32> = VecDeque::with_capacity(capacity);
    let mut vec_forward_tails: VecDeque<u32> = VecDeque::with_capacity(capacity);
    let mut vec_reverse_tails: VecDeque<u32> = VecDeque::with_capacity(capacity);
    let mut vec_mq0:           VecDeque<u32> = VecDeque::with_capacity(capacity);
    let mut vec_softclips:     VecDeque<u32> = VecDeque::with_capacity(capacity);
    let mut vec_fwd_depth:     VecDeque<u32> = VecDeque::with_capacity(capacity);
    let mut vec_rev_depth:     VecDeque<u32> = VecDeque::with_capacity(capacity);

    let mut vec_keep: VecDeque<u32> = VecDeque::with_capacity(capacity);

    match tname.find(":") {
        Some(p) => {
            bam.fetch(tname.as_str());
            tname.split_off(p);
        },
        None => {
            match header.tid(tname.as_bytes()) {
                Some(tid) => {
                    bam.fetch((tname.as_str(), 0, header.target_len(tid).unwrap()));
                },
                None => {
                    eprintln!("invalid chromosome {}", tname);
                    process::exit(0x01);
                }
            }
        }
    }
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let mut forward_heads = 0;
        let mut forward_tails = 0;
        let mut reverse_heads = 0;
        let mut reverse_tails = 0;
        let mut rev_depth = 0;
        let mut fwd_depth = 0;
        let mut softclips = 0;
        let mut mq0 = 0;
        let mut secondary = 0;
        for alignment in pileup.alignments() {
            if alignment.record().mapq() == 0 {
                mq0 += 1
            }
            if alignment.record().is_secondary() || alignment.record().mapq() == 0 {
                secondary += 1;
                continue;
            }
            // mq filter
            if alignment.record().mapq() < mqfilter {
                continue;
            }
            if alignment.record().is_reverse() {
                rev_depth += 1;
            }else {
                fwd_depth += 1;
            }
            let cigar = alignment.record().cigar();
            if alignment.is_head() {
                if cigar.leading_hardclips() > 0 || cigar.leading_softclips() > 3 || leading_insertions(&alignment.record()) > 0 {
                    softclips = softclips + 1;
                    continue;
                }
                match alignment.record().is_reverse() {
                    true => reverse_heads += 1,
                    false => forward_heads += 1,
                }
            }else if alignment.is_tail() {
                if cigar.trailing_hardclips() > 0 || cigar.trailing_softclips() > 3 || trailing_insertions(&alignment.record()) > 0 {
                    softclips = softclips + 1;
                    continue;
                }
                match alignment.record().is_reverse() {
                    true => reverse_tails += 1,
                    false => forward_tails += 1,
                }
            }
        }

        vec_depth.push_back(pileup.depth());
        vec_pos.push_back(pileup.pos());
        vec_forward_heads.push_back(forward_heads);
        vec_forward_tails.push_back(forward_tails);
        vec_reverse_heads.push_back(reverse_heads);
        vec_reverse_tails.push_back(reverse_tails);
        vec_mq0.push_back(mq0);
        vec_fwd_depth.push_back(fwd_depth);
        vec_rev_depth.push_back(rev_depth);
        vec_softclips.push_back(softclips);

        if vec_depth.len() > capacity {
            vec_depth.pop_front();
            vec_pos.pop_front();
            vec_forward_heads.pop_front();
            vec_forward_tails.pop_front();
            vec_reverse_heads.pop_front();
            vec_reverse_tails.pop_front();
            vec_mq0.pop_front();
            vec_softclips.pop_front();
            vec_fwd_depth.pop_front();
            vec_rev_depth.pop_front();
        }
        // chr22   32563973        111     20      0       0       2       0       89

        if (reverse_tails > 3 || forward_heads > 3) && vec_depth.len() > 0 {
            // skip anomalies
            if reverse_tails > 3 && (reverse_heads > reverse_tails || reverse_tails < forward_tails) {
                // eprintln!("skip at type 1 anomaries at {}", pileup.pos());
                continue;
            }else if forward_heads > 3 && (forward_heads < reverse_heads || forward_heads < forward_tails)  { // large forward heads
                // eprintln!("skip at type 2 anomaries at {}", pileup.pos());
                continue;
            }
            // push hit
            vec_keep.push_back(pileup.pos());
        }

        if vec_keep.len() > 0 && vec_keep.front().unwrap() <= &(pileup.pos() - read_length - margin) {
            let pos = vec_keep.pop_front();
            match pos {
                Some(p) => {
                    println!("---- len={}, at {}", vec_depth.len(), p);
                    for x in 0..vec_depth.len() {
                        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            tname, vec_pos[x],  vec_depth[x], vec_forward_heads[x], vec_forward_tails[x], vec_reverse_heads[x], vec_reverse_tails[x], vec_mq0[x], vec_softclips[x], vec_fwd_depth[x], vec_rev_depth[x]);
                    }
                },
                None => (println!("bad")),
            }
        }
        println!("//");
        stdout().flush();
    }
    process::exit(0x00);
}

fn leading_insertions(rec: &Record) -> i64 {
    rec.cigar().first().map_or(0, |cigar| {
        if let Cigar::Ins(s) = cigar {
            *s as i64
        }else {
            0
        }
    })
}
fn trailing_insertions(rec: &Record) -> i64 {
    rec.cigar().first().map_or(0, |cigar| {
        if let Cigar::Ins(s) = cigar {
            *s as i64
        }else {
            0
        }
    })
}
