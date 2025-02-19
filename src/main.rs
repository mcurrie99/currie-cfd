use currie_cfd as cfd;
// use cfd::matrix::Matrix;

// Constants for calculations
const IL: usize = 100001;
const N: usize = 100001;

#[allow(unused_variables)]
fn main() {
    // Creates tspan creations
    let tstart = 0.0;
    let tend = 1.5;
    let tspan = cfd::tools::linspace(tstart, tend, N);
    let dt = (tend - tstart) / ((N as f64) - 1.0);

    // Creates xspan for problems
    let xstart = 0.0;
    let xend = 5.0;
    let xspan = cfd::tools::linspace(xstart, xend, IL);
    let dx = (xend - xstart) / ((IL as f64) - 1.0);

    // Checks entered parameters will create a good soution
    let check = 2.0 * dt / dx;
    println!("Selected dt: {}", dt);
    println!("Selected dx: {}", dx);
    println!("2*dt/dx: {}", check);
    if check > 1.0 {
        panic!("ERROR: 2*dt/dx is not less than 1.0");
    }

    // Allocates memory for solutions
    let mut sol2: Vec<Vec<f64>> = vec![vec![0.0; IL]; N];
    let mut sol3: Vec<Vec<f64>> = vec![vec![0.0; IL]; N];

    // Initial Guess setup Problem 2
    for (i, val) in xspan.iter().enumerate() {
        // Sets value of the initial guess set up by equation
        if *val < 1.0 {
            sol2[0][i] = 2.0 - val;
        } else {
            sol2[0][i] = 1.0;
        }
    }

    // Initial Guess setup Problem 3
    for (i, val) in xspan.iter().enumerate() {
        // Sets value of the initial guess set up by equation
        if *val < 1.0 {
            sol3[0][i] = 2.0 - val;
        } else {
            sol3[0][i] = 1.0;
        }
    }

    // Sets up Boundary Conditions for Problem 2
    for val in sol2.iter_mut() {
        val[0] = 2.0;
    }

    // Sets up Boundary Conditions for Problem 3
    for val in sol3.iter_mut() {
        val[0] = 2.0;
    }

    // Solver Problems PDE's
    solve_pde_basic(&prob2, &mut sol2, dt, dx);
    solve_pde_basic(&prob3, &mut sol3, dt, dx);

    // Prepares to export values
    println!("Integration Complete, Exporting Results...");

    // Exporting Values for Time Discretization
    println!("Exporitng Time Grid...");
    let output_tspan = format!("tspan_{IL}_{check}.csv");
    let cols_tspan = [String::from("t")];
    let _check = cfd::tools::array_to_csv(&[tspan.clone()], &cols_tspan, &output_tspan);

    // Exports Values for Space Discretization
    println!("Exporting Space Grid...");
    let output_xspan = format!("xspan_{IL}_{check}.csv");
    let cols_xspan = [String::from("x")];
    let _check = cfd::tools::array_to_csv(&[xspan.clone()], &cols_xspan, &output_xspan);

    // Exports Calculated Values for Problem 2
    println!("Exporting Solution: Problem 2...");
    let output_prob2 = format!("problem2_{IL}_{check}.csv");
    let mut cols_prob2: Vec<String> = Vec::with_capacity(tspan.len());
    for i in 0..tspan.len() {
        let header = format!("t{}", i);
        cols_prob2.push(header);
    }
    let _check = cfd::tools::array_to_csv(&sol2, &cols_prob2, &output_prob2);

    // Exports Calculated valeus for Problem 3
    println!("Exporting Solution: Problem 3...");
    let output_prob3 = format!("problem3_{IL}_{check}.csv");
    let mut cols_prob3: Vec<String> = Vec::with_capacity(tspan.len());
    for i in 0..tspan.len() {
        let header = format!("t{}", i);
        cols_prob3.push(header);
    }
    let _check = cfd::tools::array_to_csv(&sol3, &cols_prob3, &output_prob3);

    // Outputs
    // println!("{:?}", tspan);
    // println!("{:?}", xspan);
    // println!("{:?}", sol2);
    // println!("{:?}", sol3);

}

#[allow(dead_code)]
fn solve_pde_basic<F>(prob:&F, data:&mut Vec<Vec<f64>>, dt:f64, dx:f64)
    where F: Fn(&mut Vec<Vec<f64>>, usize, usize, f64, f64) {
    // Solve Partial Differential Equation when simple
    for n in 0..data.len()-1 {
        for i in 1..data[0].len() {
            prob(data, n, i, dx, dt)
        }
    }
}

#[allow(dead_code, unused_variables)]
fn prob2(mat:&mut Vec<Vec<f64>>, ti:usize, xi:usize, dx:f64, dt:f64) {
    // First Index is Time
    // Second Index is Space
    mat[ti+1][xi] = - mat[ti][xi] * (mat[ti][xi] - mat[ti][xi-1]) / dx * dt + mat[ti][xi];
}

#[allow(dead_code, unused_variables)]
fn prob3(mat:&mut Vec<Vec<f64>>, ti:usize, xi:usize, dx:f64, dt:f64) {
    // First Index is TIme
    // Second Index is Space
    mat[ti+1][xi] = - (mat[ti][xi].powf(2.0) - mat[ti][xi-1].powf(2.0)) / (2.0 * dx) * dt + mat[ti][xi];
}
