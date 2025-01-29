// This function works with solving only single values
// NOTE: This function also approximates the current derivative by iterating small numbers 
// check for any change

// TODO: Finsih solver
#[allow(unused_variables, dead_code)]
pub fn newton_raph<F>(f:&F, x0:f64, steps:usize) -> f64 
where 
    F: Fn(f64) -> f64
{
    let mut sol = x0;

    // Iterator Solver
    for _ in 0..steps {
        // Generates new solution each step
        sol = sol - f(sol) / estimate_deriv(f, sol);
        println!("{sol}");
    }

    sol
}

// TODO: Check the functionality of this function
#[allow(unused_variables, dead_code)]
fn estimate_deriv<F>(f:&F, x:f64) -> f64 
where 
    F: Fn(f64) -> f64
{
    // Small Change in function and then estimates the slope in that region
    let change = 0.0001;
    let norm = f(x);
    let upper = f(x+change);
    let lower = f(x-change);
    let deriv_upper = (upper - norm) / change;
    let deriv_lower = (norm - lower) / change;
    let deriv=  0.5 * (deriv_upper + deriv_lower);
    deriv
}


// TODO: Make tests for these functions
#[cfg(test)]
#[allow(dead_code)]
mod test {
    use super::*;

    fn estimate_deriv_examp(x:f64) -> f64{
        let sol = 2.0*x;
        sol
    }

    #[test]
    fn test_estimate_deriv() {
        // Processes function
        let test_x = 5.0;
        let comp = estimate_deriv(&estimate_deriv_examp, test_x);
        let right = 2.0;

        // Checks if the function is working correctly
        let tolerance = 0.00001;
        let diff = (comp - right).abs();
        let check = diff < tolerance;
        assert!(check);
    }

    fn newton_raph_examp(x:f64) -> f64 {
        let sol = 3.0 - 3.0 / 40.0 * x.powf(2.0);
        sol
    }

    #[test]
    fn test_newton_raph() {
        // Test Parameters
        let tolerance = 0.0001;
        let right = 6.5;
        let steps = 1;
        let init_guess = 5.0;

        // Runs the test and compares
        let test = newton_raph(&newton_raph_examp, init_guess, steps);
        let diff = test - right;
        let check = diff < tolerance;
        assert!(check);
    }
}
