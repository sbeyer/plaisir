use std::env::current_dir;
use std::process::Command;

fn main() {
    let cwd = current_dir().expect("get current working directory failed");
    Command::new("make")
        .args(&["-C", "lkh"])
        .status()
        .expect("make failed");
    println!("cargo:rustc-link-search={}/lkh/", cwd.display());
    println!("cargo:rustc-link-lib=LKH");
}
