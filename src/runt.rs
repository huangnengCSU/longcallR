#[derive(Default, Debug, Clone)]
pub struct Runtime{
    pub readname: String,
    pub length: u32,
    pub runtime: u64,
}

impl Runtime {
    pub fn new(readname: String, length: u32, runtime: u64) -> Runtime {
        Runtime {
            readname,
            length,
            runtime,
        }
    }
}