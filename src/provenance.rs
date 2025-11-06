use std::fs::File;
use std::hash::Hasher;
use std::io::{Read, Result as IoResult};
use std::path::Path;
use twox_hash::XxHash64;

pub fn filehash_xx64(path: &Path) -> IoResult<String> {
    let mut file = File::open(path)?;
    let mut hasher = XxHash64::with_seed(0);
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.write(&buf[..n]);
    }
    Ok(format!("{:016x}", hasher.finish()))
}
