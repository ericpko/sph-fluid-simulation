use std::error::Error;


fn main() -> Result<(), Box<dyn Error>> {
   sph_fluid_simulation::simulate()?;
   Ok(())
}
