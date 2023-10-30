#[derive(Debug, Default, Clone)]
pub struct VCFRecord {
    pub chromosome: Vec<u8>,
    pub position: u64,
    pub id: Vec<u8>,
    pub reference: Vec<u8>,
    pub alternative: Vec<Vec<u8>>,
    pub qual: f32,
    pub filter: Vec<u8>,
    pub info: Vec<u8>,
    pub format: Vec<u8>,
    pub genotype: String,
    /* private fields */
}