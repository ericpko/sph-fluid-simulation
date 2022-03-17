use glam::Vec2;
use rand::Rng;

#[derive(Clone)]
pub struct Particle {
    pub pos: Vec2,
    pub vel: Vec2,
    pub force: Vec2,
    pub rho: f32,
    pub pressure: f32,
}

impl Particle {
    pub fn new(win_width: f32, win_height: f32) -> Self {
        let mut rng = rand::thread_rng();
        let pos = Vec2::new(
            rng.gen_range(
                (win_width / 2.0 - (win_width / 4.0))..(win_width / 2.0 + (win_width / 4.0)),
            ),
            rng.gen_range(
                (win_height / 2.0 - (win_height / 4.0))..(win_height / 2.0 + (win_height / 4.0)),
            ),
        );
        let vel = Vec2::new(0.0, 0.0);
        let force = Vec2::new(0.0, 0.0);

        Self {
            pos,
            vel,
            force,
            rho: 0.0,
            pressure: 0.0,
        }
    }
}
