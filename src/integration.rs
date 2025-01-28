// TODO: Implement Structure to have easy retreival and export of integration data
struct Integration {
    t: Vec<f64>,
    y: Vec<Vec<f64>>,
    int_type: IntegrationType,

}

// TODO: Needs more fullfiled implementation
// To inform user of the type of integration that was used
enum IntegrationType {
    RK2,
    AM2,
}

// 64 bit module for numerical integration
// RUNGE-KUTTA 2ND ORDER METHOD
pub fn rk2<F>(f: &F, x0:&Vec<f64>, tspan:&Vec<f64>) -> Vec<Vec<f64>> 
where
    F: Fn(&Vec<f64>, f64) -> Vec<f64>
{
    // f: Function to call for numerical integration
    // x0: Initial Conditions
    // tspan: Integration Axis
    let mut results= vec![vec![0.0; tspan.len()]; x0.len() + 1];
    
    // Initial Calculations using x0
    // Adds time to calculations
    // Iterates through initial variables before beginning integration
    results[0][0] = *tspan.get(0).unwrap();
    for (i, val) in x0.iter().enumerate() {
        results[i+1][0] = *val;
    }
    
    // Allocates input variables
    let mut inputs = vec![0.0; x0.len()];

    // Analysis over time span
    for (i, time) in tspan.iter().enumerate().take(tspan.len() - 1) {
        
        // Constructs inputs
        for j in 1..results.len() {
            inputs[j-1] = results[j][i];
        }
        
        // First Stage of Analysis
        let dt = tspan.get(i+1).unwrap() - tspan.get(i).unwrap();
        let (guess, derivs) = rk2_guess(&f, &inputs, *time, dt);

        // Second Stage of analysis
        let vals = rk2_act(&f, &guess, *time, dt, &derivs, &inputs);

        // Enters values into the results array
        results[0][i+1] = tspan[i+1];
        for (j, val) in vals.iter().enumerate() {
            results[j+1][i+1] = *val;
        }
    }

    results
    
}

// Initial Guess for the integration
fn rk2_guess<F>(f: &F, inputs:&Vec<f64>, time:f64, dt:f64) -> (Vec<f64>, Vec<f64>) 
where 
    F: Fn(&Vec<f64>, f64) -> Vec<f64>
{
    // Allocates memory for results
    let mut results = vec![0.0; inputs.len()];
    
    // Call Derivative funciton
    let derivs = f(&inputs, time);

    // Ensures that the bounds of both arrays are equal
    assert_eq!(&inputs.len(), &derivs.len());

    // Goes through and asserts what the guesses are
    for (i, deriv) in derivs.iter().enumerate() {
        results[i] = deriv * dt + inputs[i];
    }

    (results, derivs)
}

fn rk2_act<F>(f: &F, guess:&Vec<f64>, time:f64, dt:f64, derivs:&Vec<f64>, inputs:&Vec<f64>) -> Vec<f64> 
where 
    F: Fn(&Vec<f64>, f64) -> Vec<f64>
{
    // Allocates the appropriate memory
    let mut results = vec![0.0; guess.len()];
    
    // Calls derivative on current values
    let derivs_new = f(&guess, time+dt);

    // Asserts that the lenghts are the same
    assert_eq!(&derivs_new.len(), &guess.len());

    // Calculates values
    for (i, deriv_new) in derivs_new.iter().enumerate() {
        results[i] = 0.5 * (derivs[i] + deriv_new) * dt + inputs.get(i).unwrap();
    }

    results
}


// TODO: Finish this, I have grown tired of my state
// I will defeat this one day
// 2ND ORDER ADAMS-MOULTON
#[allow(unused_variables)]
pub fn am2<F>(f:&F, x0:&Vec<f64>, tspan:&Vec<f64>, steps:usize) -> Vec<Vec<f64>> 
where F:Fn(&Vec<f64>, f64) -> Vec<f64>{
    // f: Function to call for numerical integration
    // x0: Initial Conditions
    // tspan: Integration Axis
    // steps: Number of Steps the linear Solver should make
    println!("{}", x0.len());
    let mut results= vec![vec![0.0; tspan.len()]; x0.len() + 1];
    println!("{}", results.len());

    // Initial Calculations using x0
    // Adds time to calculations
    // Iterates through initial variables before beginning integration
    results[0][0] = *tspan.get(0).unwrap();
    for (i, val) in x0.iter().enumerate() {
        results[i+1][0] = *val;
    }

    // Allocates input variables
    let mut inputs = vec![0.0; x0.len()];

    // Sets results equal to a clone of tspan
    results[0] = tspan.clone();

    // Analysis over time span
    for (i, time) in tspan.iter().enumerate().take(tspan.len() - 1) {        
        // Constructs inputs
        for j in 1..results.len() {
            inputs[j-1] = results[j][i];
        }

        println!("{:?}", inputs);

        // Sets up some variables
        let dt = tspan.get(i+1).unwrap() - tspan.get(i).unwrap();

        // Creates Closure of root function and solves
        // Enters values into the results array
        for j in 0..x0.len() {
            println!("{}", j);
            // Solves values for the timestep
            let solve = |yn1| 
                0.5 * (f(&inputs, *time).get(j).unwrap() + f(&vec![yn1], *time).get(j).unwrap()) 
                - (yn1 - inputs.get(j).unwrap()) / dt;
            
            // Solves presented ODE
            let ans = super::linear::newton_raph(&solve, 5.0, steps);
            results[j+1][i+1] = ans;

        }
    }

    results
}

// TODO: I think this can be removed
#[allow(unused_variables, dead_code)]
fn am2_roots(yn0:f64, yn1:f64) -> f64 {
    5.0
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper Function for Asserting tolerance
    fn asssert_eq_tol(a:f64, b:f64, tolerance:f64) {
        let difference = (a - b).abs();
        assert!(
            difference <= tolerance,
            "Integration Failed to meet acceptable tolerance"
        );
    }

    // Derivative Function used for test
    fn derivs_int(x: &Vec<f64>, t: f64) -> Vec<f64> {
        // Assess current values
        let y = x.get(0).unwrap();
        let z = x.get(1).unwrap();
        let x = t;
    
        let dy_dx = z.clone();
    
        let dy2_dx2 = 5.0 * x.powf(2.0) - 3.0 / 2.0 * x.powf(2.0) * y * z;
        
        vec![dy_dx, dy2_dx2]
    }   

    #[test]
    fn test_rk2() {
        // Defines Tolerance for Test Case
        let tolerance = 0.0001;

        // Defines Solution to Test Case
        let sol = vec![
            vec![1.0, 1.1, 1.2],
            vec![2.0, 2.195, 2.37525],
            vec![2.0, 1.873165, 1.716934088]
        ];

        // Sets up Initial Conditions for test
        let y0 = 2.0;
        let z0 = 2.0;

        // Inputs for integration
        let input = vec![y0, z0];
        let tspan = vec![1.0, 1.1, 1.2];

        // Calculates Results from Integration
        let int_results = rk2(&derivs_int, &input, &tspan);


        for (row_sol, row_calc) in sol.iter().zip(int_results.iter()) {
            for (comp_sol, comp_calc) in row_sol.iter().zip(row_calc.iter()) {
                asssert_eq_tol(*comp_sol, *comp_calc, tolerance);
            }
        }
    }
}