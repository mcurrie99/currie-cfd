// Crate that is being worked on
use currie_cfd as cfd;

#[allow(dead_code)]
fn main2() {
    let y0 = 2.0;
    let z0 = 2.0;
    let x0 = 1.0;

    // let test = Float::with_val(128, 0.1+0.2);
    // println!("{}", test);

    let input = vec![y0, z0];
    let tinput = x0;
    
    // Creates tspan Vector
    let dt = 0.0001;
    let mut tspan: Vec<f64> = Vec::with_capacity(1000); 
    tspan.push(1.0);
    for i in 1..2000 {
        tspan.push(tspan[i-1]+dt);
    }
    
    // let tspan = vec![1.0, 1.1, 1.2];

    let result = test_int(&input, tinput);
    println!("{:?}", result);

    let int_results = cfd::integration::rk2(&test_int, &input, &tspan);

    println!("{:?}", int_results);
    
    let output_path = "output2.csv";
    let cols = [
        String::from("Time"), 
        String::from("Y"), 
        String::from("Z")];
    let _check = cfd::tools::array_to_csv(&int_results, &cols, output_path);
}

fn main() {
    // Setting up conditions
    let y0 = 10.0;
    
    // Sets up Arrays for funciton call
    let xspan = vec![1.0, 1.1, 1.2];
    let input = vec![y0];

    // Defining number of steps for solver
    let solver_steps = 100;

    // Creates tspan Vector
    // let dt = 0.0001;
    // let mut xspan: Vec<f64> = Vec::with_capacity(1000); 
    // xspan.push(1.0);
    // for i in 1..2000 {
    //     xspan.push(xspan[i-1]+dt);
    // }

    // Calculates the results
    let int_results = cfd::integration::am2(&test_again, &input, &xspan, solver_steps);

    // Prints results to the console
    println!("{:?}", int_results);

    let output_path = "am2_act.csv";
    let cols = [
        String::from("Time"), 
        String::from("Y")];
    let _check = cfd::tools::array_to_csv(&int_results, &cols, output_path);
}

#[allow(dead_code)]
fn test_int(x: &Vec<f64>, t: f64) -> Vec<f64> {
    // Assess current values
    let y = x.get(0).unwrap();
    let z = x.get(1).unwrap();
    let x = t;

    let dy_dx = z.clone();

    let dy2_dx2 = 5.0 * x.powf(2.0) - 3.0 / 2.0 * x.powf(2.0) * y * z;
    
    vec![dy_dx, dy2_dx2]
}   

// TODO: Remove this when you have completed the homework
#[allow(unused_variables, dead_code)]
fn test_again(x: &Vec<f64>, t:f64) -> Vec<f64> {
    // Assess current values
    let yn = x.get(0).unwrap();

    let dy_dx = 5.0 - 3.0 / 2.0 * yn.powf(2.0);

    vec![dy_dx]
}

