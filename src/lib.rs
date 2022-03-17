/**
 * Particle-Based Fluid Simulation for Interactive Applications
 * https://matthias-research.github.io/pages/
 */
use ggez::event;
use ggez::graphics::{self, Color, Mesh};
use ggez::{Context, GameResult};
use glam::Vec2;
use std::f32::consts::PI;

mod particle;
use particle::Particle;

// Window dimensions
// const ASPECT_RATIO: f32 = 16.0 / 9.0;
// const WIDTH: f32 = 1280.0;
// const HEIGHT: f32 = WIDTH / ASPECT_RATIO;
const WIDTH: f32 = 350.0;
const HEIGHT: f32 = 600.0;
const WIN_PADDING: f32 = 20.0;

// Simulation constants
const N_PARTICLES: usize = usize::pow(20, 2);
const MASS: f32 = 10.0;
const _DT: f32 = 1e-3;
const H: f32 = 2.8; // Radius of support
const GRAVITY: f32 = -1e4; // Gravity
const GAS_K: f32 = 1350.0; // Gas constant k
const RHO_0: f32 = 1e3; // Rest density for water
const MU: f32 = 160.0; // Viscosity
const DAMPENING: f32 = 0.8; // Wall collision dampening

// Smoothing kernels and derivatives
const POLY6: f32 = 315.0 / (64.0 * PI * H * H * H * H * H * H * H * H * H);
const D_SPIKY: f32 = -45.0 / (PI * H * H * H * H * H * H);
const L_VISC: f32 = 45.0 / (PI * H * H * H * H * H * H);

enum PlayState {
    Play,
    Pause,
    Restart
}

struct ParticleSystemState {
    state: PlayState,
    mesh: Mesh,
    dt: std::time::Duration,
    particles: Vec<Particle>,
    neighbors: Vec<Vec<usize>>,
}

impl ParticleSystemState {
    fn new(ctx: &mut Context) -> GameResult<ParticleSystemState> {
        let state = PlayState::Pause;
        let mesh = graphics::Mesh::new_circle(
            ctx,
            graphics::DrawMode::fill(),
            Vec2::new(0.0, 0.0),
            5.0,
            0.3,
            Color::from_rgb(0, 255, 255),
        )?;
        let dt = std::time::Duration::new(0, 0);
        let particles = Self::init_particles();
        let neighbors = Self::init_neighbors();
        let s = ParticleSystemState {
            state,
            mesh,
            dt,
            particles,
            neighbors,
        };
        Ok(s)
    }

    fn init_particles() -> Vec<Particle> {
        std::iter::repeat_with(|| Particle::new(WIDTH, HEIGHT))
            .take(N_PARTICLES)
            .collect()
    }

    fn init_neighbors() -> Vec<Vec<usize>> {
        let mut neighbors = Vec::with_capacity(N_PARTICLES);
        for _ in 0..N_PARTICLES {
            let neighbors_i: Vec<usize> = Vec::new();
            neighbors.push(neighbors_i);
        }
        return neighbors;
    }

    fn reset_particles(&mut self) {
        self.particles = Self::init_particles();
    }

    fn update_neighbors(&mut self) {
        let h2 = H * H;
        for i in 0..self.particles.len() {
            // let mut neighbors_i = Vec::new();
            self.neighbors[i].clear();
            let r_i = self.particles[i].pos;
            for j in 0..self.particles.len() {
                let r_j = self.particles[j].pos;
                let dx = r_j.x - r_i.x;
                let dy = r_j.y - r_i.y;
                let r2 = dx * dx + dy * dy;
                if r2 < h2 {
                    self.neighbors[i].push(j);
                }
            }
        }
    }

    fn zero_density_pressure_force(&mut self) {
        for p in &mut self.particles {
            p.force = Vec2::new(0.0, 0.0);
            p.rho = 0.0;
            p.pressure = 0.0;
        }
    }

    fn compute_density_and_pressure(&mut self) {
        let h2 = H * H;
        for i in 0..self.particles.len() {
            let r_i = self.particles[i].pos;
            for j in &self.neighbors[i] {
                let r_j = self.particles[*j].pos;
                let dx = r_j.x - r_i.x;
                let dy = r_j.y - r_i.y;
                let r2 = dx * dx + dy * dy;
                self.particles[i].rho += MASS * POLY6 * f32::powf(h2 - r2, 3.0);
            }
            self.particles[i].pressure = GAS_K * (self.particles[i].rho - RHO_0);
        }
    }

    fn compute_forces(&mut self) {
        for i in 0..self.particles.len() {
            self.particles[i].force.y = -GRAVITY;
            self.neighbors[i].retain(|&x| x != i);
            if self.neighbors[i].len() == 0 {
                continue;
            }
            let r_i = self.particles[i].pos;
            let u_i = self.particles[i].vel;
            let mut force_pressure_i = Vec2::new(0.0, 0.0);
            let mut force_visc_i = Vec2::new(0.0, 0.0);
            for j in &self.neighbors[i] {
                let r_j = self.particles[*j].pos;
                let u_j = self.particles[*j].vel;
                let dx = r_j.x - r_i.x;
                let dy = r_j.y - r_i.y;

                let r = f32::sqrt(dx * dx + dy * dy);
                let rhat = (r_j - r_i) / r;

                force_pressure_i += -MASS
                    * ((self.particles[i].pressure + self.particles[*j].pressure)
                        / (2.0 * self.particles[*j].rho))
                    * D_SPIKY
                    * f32::powf(H - r, 2.0)
                    * rhat;
                force_visc_i +=
                    MU * MASS * ((u_j - u_i) / self.particles[*j].rho) * L_VISC * (H - r);
            }
            self.particles[i].force += force_pressure_i + force_visc_i;
        }
    }

    fn symplectic_euler(&mut self, _dt: f32) {
        let left = 0.0 + WIN_PADDING;
        let right = WIDTH - WIN_PADDING;
        let bottom = 0.0 + WIN_PADDING;
        let top = HEIGHT - WIN_PADDING;

        for i in 0..self.particles.len() {
            let mut u_i =
                self.particles[i].vel + _DT * self.particles[i].force / self.particles[i].rho;
            let mut r_i = self.particles[i].pos + _DT * u_i;

            if r_i.x < left || r_i.x > right {
                u_i.x = -u_i.x * DAMPENING;
                r_i.x = if r_i.x < left { left } else { right }
            }

            if r_i.y < bottom || r_i.y > top {
                u_i.y = -u_i.y * DAMPENING;
                r_i.y = if r_i.y < bottom { bottom } else { top }
            }
            self.particles[i].pos = r_i;
            self.particles[i].vel = u_i;
        }
    }
}

impl event::EventHandler<ggez::GameError> for ParticleSystemState {
    fn update(&mut self, ctx: &mut Context) -> GameResult {
        self.dt = ggez::timer::delta(ctx);
        let tick = self.dt.as_secs_f32();

        let pressed_key = ggez::input::keyboard::pressed_keys(ctx);
        if pressed_key.contains(&ggez::event::KeyCode::Space) {
            self.state = PlayState::Play;
        } else if pressed_key.contains(&ggez::event::KeyCode::R) {
            self.state = PlayState::Restart;
        } else if pressed_key.contains(&ggez::event::KeyCode::P) {
            self.state = PlayState::Pause;
        }

        match self.state {
            PlayState::Play => {
                self.update_neighbors();
                self.zero_density_pressure_force();
                self.compute_density_and_pressure();
                self.compute_forces();
                self.symplectic_euler(tick);
            }
            PlayState::Pause => {
            }
            PlayState::Restart => {
                self.reset_particles();
                self.state = PlayState::Pause;
            }
        }

        Ok(())
    }

    fn draw(&mut self, ctx: &mut Context) -> GameResult {
        graphics::clear(ctx, Color::from_rgb(30, 30, 30));

        // let circle = graphics::Mesh::new_circle(
        //     ctx,
        //     graphics::DrawMode::fill(),
        //     Vec2::new(0.0, 0.0),
        //     5.0,
        //     0.3,
        //     Color::from_rgb(0, 255, 255)
        // )?;
        for p in &self.particles {
            let params = graphics::DrawParam::new().dest(p.pos);
            // graphics::draw(ctx, &circle, params)?;
            graphics::draw(ctx, &self.mesh, params)?;
        }

        graphics::present(ctx)?;
        Ok(())
    }
}

pub fn simulate() -> GameResult {
    let mut cb = ggez::ContextBuilder::new("SPH-Fluid", "Eric Koehli");

    if let Ok(manifest_dir) = std::env::var("CARGO_MANIFEST_DIR") {
        let mut path = std::path::PathBuf::from(manifest_dir);
        path.push("resources");
        println!("Adding path {:?}", path);
        cb = cb.add_resource_path(path);
    }

    let (mut ctx, event_loop) = cb
        .window_mode(
            ggez::conf::WindowMode::default()
                .dimensions(WIDTH, HEIGHT)
                .resizable(true),
        )
        .window_setup(
            ggez::conf::WindowSetup::default()
                .title("SPH Fluid Simulation!")
                .samples(ggez::conf::NumSamples::Eight)
                .vsync(true),
        )
        .build()?;
    let state = ParticleSystemState::new(&mut ctx)?;
    println!("\n=====================\nPress <Space> to run\nPress <P> to pause\nPress <R> to reset");
    event::run(ctx, event_loop, state)
}
