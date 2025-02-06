#[allow(dead_code)]
pub struct MatrixSparse {
    pub vals: Vec<(f64, f64, f64)>,
    rows: usize,
    cols: usize
}

#[allow(dead_code)]
pub struct Matrix {
    pub vals: Vec<Vec<f64>>,
    rows: usize,
    cols: usize,
}

#[allow(dead_code)]
impl Matrix {
    pub fn new(rows: usize, cols: usize) -> Matrix {
        Matrix {
            vals: vec![vec![0.0; cols]; rows],
            rows,
            cols
        }
    }

    pub fn rows(&self) -> usize {
        self.rows.clone()
    }

    pub fn cols(&self) -> usize {
        self.cols.clone()
    }
}

// impl Add for Matrix {
//     fn add(self, Rhs:Matrix) -> Matrix {
//         // Asserts that the dimensions are equal to each other
//         assert_eq!(self.cols, Rhs.cols);
//         assert_eq!(self.rows, Rhs.rows);
    
//         // Creates new Matrix
//         let mut result = Matrix::new(self.rows, self.cols);

//         // Adds
//         for i in 0..self.rows {
//             for j in 0..self.cols {
//                 result.vals[i][j] = self.vals[i][j] + Rhs.vals[i][j];
//             }
//         }

//         result
//     }
// }


// TODO: Testing for Matrix Stuff
        // // println!("{}, {}, {}", guess[0], guess[1], guess[2]);
        // // println!("{dx}");
        // // println!("{}", mat.vals[0][0]);
        // println!("{:?}", result.vals);

        // let mut test = Matrix::new(3,3);
        // test.vals[0] = vec![2.0, 7.0, 0.0];
        // test.vals[1] = vec![3.0, -2.0, 2.0];
        // test.vals[2] = vec![0.0, 5.0, 3.0];

        // let mut tester = Matrix::new(3, 1);
        // tester.vals = vec![vec![1.0; 3]; 3];

        // println!("A:\n{:?}", test.vals);
        // println!("b:\n{:?}", tester.vals);
        // let (test_l, test_u) = build_lu(&test);

        // // Outputs Matrices
        // println!("L:\n{:?}", test_l.vals);
        // println!("U:\n{:?}", test_u.vals);

        // // Outputs solution
        // let (test_x, test_z) = solve_lu(&test_l, &test_u, &tester);
        // println!("x:\n{:?}", test_x.vals);
