// Crate that is being worked on;
use super::matrix::Matrix;
// use super::ivp;

// Constants for the problem
const A: f64 = 10.0;
const B: f64 = 1e-4;
const D: f64 = 4.0e5;
const E: f64 = 2.0e3;
const TH: f64 = 2000.0;
const TC: f64 = 400.0;
const HH: f64 = 900.0;
const HL: f64 = 1100.0;
const L: f64 = 0.1;



#[allow(unused_variables, non_snake_case)]
fn main() {
    // Setup of model
    let mut error = f64::INFINITY;
    let tolerance = 1e-10;
    let length = 10001;

    let start_guess = 1600.0;
    let end_guess = 600.0;

    // dx derivation
    let dt = L / (length - 1) as f64;

    // Sets up initial guess
    let mut guess = vec![TH; length];
    for i in 0..length {
        guess[i] = start_guess - (start_guess - end_guess) / (L) * dt * (i as f64);
    }
    println!("{:?}", guess);

    // Sets up solving matrix
    let mut guess_mat = Matrix::new(guess.len()-2, 1);
    for j in 0..guess_mat.rows() {
        guess_mat.vals[j][0] = guess[j+1];
    }

    println!("{:?}", guess_mat.vals);

    // Creates the matrix that the math will be done with
    let mut mat = Matrix::new(guess.len()-2, guess.len()-2);
    let mut result = Matrix::new(guess.len()-2, 1);
    // let mut dT = Matrix::new(guess.len()-2, 1);


    // TODO: Move this to generic Function
    // Function Calling
    while error > tolerance {
        // Assembles Matrix
        for i in 0..mat.rows() {
            // Grabs current values
            let tin1 = *guess.get(i).unwrap();
            let ti = *guess.get(i+1).unwrap();
            let tip1 = *guess.get(i+2).unwrap();

            // if i==0 {
            //     println!("PP:{}", tin1)
            // }

            // Grabs Derivatives
            let derivs = derivs(tin1, ti, tip1, dt);

            // Determines values to place into Matrix
            if i == 0 {
                // If it is the first iteration
                mat.vals[i][i] = *derivs.get(1).unwrap();
                mat.vals[i][i+1] = *derivs.get(2).unwrap();
            } else if i == mat.rows()-1 {
                // If it is the last iteration
                mat.vals[i][i-1] = *derivs.get(0).unwrap();
                mat.vals[i][i] = *derivs.get(1).unwrap();
            } else {
                // Any other part of the matrix
                mat.vals[i][i-1] = *derivs.get(0).unwrap();
                mat.vals[i][i] = *derivs.get(1).unwrap();
                mat.vals[i][i+1] = *derivs.get(2).unwrap();
            }

            // Determines Values to place into resulting vector
            result.vals[i] = vec![-func(tin1, ti, tip1, dt)];
        }   

        // Decomposes matrix and solves
        let (l, u) = build_lu(&mat);
        let (dx, z) = solve_lu(&mat, &l, &u, &result);

        // Adjusts Temperatues
        for j in 0..result.rows() {
            guess_mat.vals[j][0] = guess_mat.vals[j][0] + dx.vals[j][0];
            guess[j+1] = guess_mat.vals[j][0]
        }

        // let first_val = |t1| HH*(TH-t1) + (A + B*guess[1]) * (-3.0*t1 + 4.0*guess[1] - guess[2]) / (2.0 * dt);
        // guess[0] = cfd::linear::newton_raph(&first_val, start_guess, 14);
        // println!("Guess: {}", guess[0]);
    
        // let idx = guess.len()-1;
        // let last_val = |til:f64| HL * (til - TC) + (A + B * til) * (3.0 * til - 4.0*guess[idx-1] + guess[idx-2]) / (2.0 * dt);
        // guess[idx] = cfd::linear::newton_raph(&last_val, end_guess, 14);

        // Calculates Error
        error = find_error(&guess_mat, &dx);
        // println!("{}", error);
        // println!("dT:\n{:?}", dx.vals);
        // println!("F:\n{:?}", result.vals);
        // println!("Matrix:\n{:?}", mat.vals);
        // error = 0.0;
    
    }

    println!("Relative Error: {}", error);


    // Moves values to array of numbers
    let mut results = vec![0.0; length];
    for i in 1..results.len()-1 {
        results[i] = guess_mat.vals[i-1][0];
    }

    let first_val = |t1| HH*(TH-t1) + (A + B*results[1]) * (-3.0*t1 + 4.0*results[1] - results[2]) / (2.0 * dt);
    results[0] = crate::linear::newton_raph(&first_val, start_guess, 14);

    let idx = results.len()-1;
    let last_val = |til:f64| HL * (til - TC) + (A + B * til) * (3.0 * til - 4.0*results[idx-1] + results[idx-2]) / (2.0 * dt);
    results[idx] = crate::linear::newton_raph(&last_val, end_guess, 14);
    // println!("{:?}", results);

    // Exports values
    let mut output = vec![vec![0.0; results.len()]; 2];
    for i in 0..output[0].len() {
        output[0][i] = (i as f64) * dt;
        output[1][i] = results[i] as f64;
    }

    let output_path = "edited.csv";
    let cols = [
        String::from("x"), 
        String::from("T")];
    let _check = crate::tools::array_to_csv(&output, &cols, output_path);
}

// Function for finding derivatives
#[allow(dead_code, non_snake_case)]
fn derivs(Tin1:f64, Ti:f64, Tip1:f64, dx:f64) -> Vec<f64> {
    let deriv_Tin1 = (Tin1 * B + A) / dx.powf(2.0);
    let deriv_Ti = (-2.0 * Ti * B - 2.0 * A) / dx.powf(2.0) + E;
    let deriv_Tip1 = (Tip1 * B + A) / dx.powf(2.0);

    // Returns Values
    vec![deriv_Tin1, deriv_Ti, deriv_Tip1]
}

// Ouput of Function
#[allow(dead_code, non_snake_case)]
fn func(Tin1:f64, Ti:f64, Tip1:f64, dx:f64) -> f64 {
    let part1 = dx.powf(2.0) * (Ti * E + D);
    let part2 = -(Ti - Tip1) * (0.5 * Ti * B + 0.5 * Tip1 * B + A);
    let part3 = -(Ti - Tin1) * (0.5 * Ti * B + 0.5 * Tin1 * B + A);
    let result = (part1 + part2 + part3) / dx.powf(2.0);
    result
}

#[allow(dead_code)]
fn diff_op_1(x1:f64, x2:f64, dx:f64)  -> f64 {
    (x1 - x2) / dx
}

#[allow(dead_code)]
fn diff_op_2(x1:f64, x2:f64, x3:f64, dx:f64) -> f64 {
    (x1 - 2.0*x2 + x3) / dx.powf(2.0)
}

#[allow(dead_code, non_snake_case)]
fn test_func(Tin1:f64, Ti:f64, Tip1:f64, dx:f64) -> f64 {
    let part1 = B * diff_op_1(Ti, Tin1, dx).powf(2.0);
    let part2 = (A + B * Ti) * diff_op_2(Tip1, Ti, Tin1, dx);
    let part3 = D + E * Ti;

    part1 + part2 + part3
}

#[allow(dead_code, non_snake_case)]
fn test_derivs(Tin1:f64, Ti:f64, Tip1:f64, dx:f64) -> Vec<f64> {
    let deriv_Tin1 = (1.0/dx.powf(2.0))*(0.5*B*Tin1-0.5*B*Tip1+B*Ti+A);
    let deriv_Ti = (1.0/dx.powf(2.0))*(-4.0*B*Ti+B*Tip1+B*Tin1-2.0*A)+E;
    let deriv_Tip1 = (1.0/dx.powf(2.0))*(0.5*B*Tip1-0.5*B*Tin1+B*Ti+A);

    vec![deriv_Tin1, deriv_Ti, deriv_Tip1]
}

// Functions for solving
#[allow(dead_code)]
fn build_lu(mat: &Matrix) -> (Matrix, Matrix) {
    // Generates L and U matrices
    let mut l = Matrix::new(mat.rows(), mat.cols());
    let mut u = Matrix::new(mat.rows(), mat.cols());
    
    // Places first entries into arrays
    l.vals[0][0] = mat.vals[0][0];
    u.vals[0][0] = 1.0;
    u.vals[0][1] = mat.vals[0][1] / l.vals[0][0];

    // Calculates elements of L and U matrixes
    for i in 1..mat.rows() {
        
        // Grabs values from array
        let ci = mat.vals[i][i-1];
        let ai = mat.vals[i][i];

        // Value Creation
        let pi = ci;
        let qin1 = u.vals[i-1][i];
        let li = ai - pi * qin1;

        // Assignment to matrix
        l.vals[i][i] = li;
        l.vals[i][i-1] = pi;
        u.vals[i][i] = 1.0;

        // Assigns q if not the last row
        if i != mat.rows()-1 {
            let bi = mat.vals[i][i+1];
            let qi = bi / li;
            u.vals[i][i+1] = qi;
        }
    }
    
    // Returns the two matrices
    (l, u)
}

#[allow(dead_code)]
fn solve_lu(mat: &Matrix, l: &Matrix, u: &Matrix, b: &Matrix) -> (Matrix, Matrix) {
    // Returns the z and x Matrix objects to the user from solving
    // l: lower triangular matrix
    // u: upper triangular matrix
    // b: vector that wants to be solved
    // Returns (x, z)
    
    // Assertions to ensure assumptions of incoming matrices
    assert_eq!(l.rows(), u.rows());
    assert_eq!(l.cols(), u.cols());
    assert_eq!(l.cols(), b.rows());

    let mut z = Matrix::new(b.rows(), 1);
    let mut x = Matrix::new(b.rows(), 1);

    // Enters first value
    z.vals[0][0] = b.vals[0][0] / l.vals[0][0];

    // Enters remaining values
    for i in 1..l.rows() {
        z.vals[i][0] = (b.vals[i][0] - mat.vals[i][i-1] * z.vals[i-1][0]) / l.vals[i][i];
    }

    // Starts Assesing the x vector

    // First entry into x vector
    let len = x.rows()-1;
    x.vals[len] = vec![z.vals[z.rows()-1][0]];

    // Following entries in to the vector
    for i in 1..=len {
        let idx = len - i;
        x.vals[idx][0] = z.vals[idx][0] - u.vals[idx][idx+1] * x.vals[idx+1][0];
    }

    // Returns vectors
    (x, z)
}


fn find_error(x:&Matrix, dx:&Matrix) -> f64 {
    // Makes assertions
    assert_eq!(x.rows(), dx.rows());

    let mut error: f64 = 0.0;

    for i in 0..x.rows() {
        // Adds to error
        let num = dx.vals[i][0] as f64;
        let den = x.vals[i][0] as f64;
        error += (num / den).powf(2.0);
    }

    // Takes square root of error
    error = error.sqrt();

    // Returns values
    error
}


// #[allow(dead_code)]
// fn main2() {
//     let y0 = 2.0;
//     let z0 = 2.0;
//     let x0 = 1.0;

//     // let test = Float::with_val(128, 0.1+0.2);
//     // println!("{}", test);

//     let input = vec![y0, z0];
//     let tinput = x0;
    
//     // Creates tspan Vector
//     let dt = 0.0001;
//     let mut tspan: Vec<f64> = Vec::with_capacity(1000); 
//     tspan.push(1.0);
//     for i in 1..2000 {
//         tspan.push(tspan[i-1]+dt);
//     }
    
//     // let tspan = vec![1.0, 1.1, 1.2];

//     let result = test_int(&input, tinput);
//     println!("{:?}", result);

//     let int_results = cfd::ivp::rk2(&test_int, &input, &tspan);

//     println!("{:?}", int_results);
    
//     let output_path = "output2.csv";
//     let cols = [
//         String::from("Time"), 
//         String::from("Y"), 
//         String::from("Z")];
//     let _check = cfd::tools::array_to_csv(&int_results, &cols, output_path);
// }

// fn main() {
//     // Setting up conditions
//     let y0 = 10.0;
    
//     // Sets up Arrays for funciton call
//     let xspan = vec![1.0, 1.1, 1.2];
//     let input = vec![y0];

//     // Defining number of steps for solver
//     let solver_steps = 2;

//     // Creates tspan Vector
//     // let dt = 0.0001;
//     // let mut xspan: Vec<f64> = Vec::with_capacity(1000); 
//     // xspan.push(1.0);
//     // for i in 1..2000 {
//     //     xspan.push(xspan[i-1]+dt);
//     // }

//     // Calculates the results
//     let int_results = cfd::ivp::am2(&test_again, &input, &xspan, solver_steps);

//     // Prints results to the console
//     println!("{:?}", int_results);

//     let output_path = "am2_act.csv";
//     let cols = [
//         String::from("Time"), 
//         String::from("Y")];
//     let _check = cfd::tools::array_to_csv(&int_results, &cols, output_path);
// }

// #[allow(dead_code)]
// fn test_int(x: &Vec<f64>, t: f64) -> Vec<f64> {
//     // Assess current values
//     let y = x.get(0).unwrap();
//     let z = x.get(1).unwrap();
//     let x = t;

//     let dy_dx = z.clone();

//     let dy2_dx2 = 5.0 * x.powf(2.0) - 3.0 / 2.0 * x.powf(2.0) * y * z;
    
//     vec![dy_dx, dy2_dx2]
// }

// // TODO: Remove this when you have completed the homework
// #[allow(unused_variables, dead_code)]
// fn test_again(x: &Vec<f64>, t:f64) -> Vec<f64> {
//     // Assess current values
//     let yn = x.get(0).unwrap();

//     let dy_dx = 5.0 - 3.0 / 2.0 * yn.powf(2.0);

//     vec![dy_dx]
// }

