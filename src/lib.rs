// Integration Module
pub mod ivp;

// Linear Solver Module
pub mod linear;

// Matrix Implementation
pub mod matrix;

// Boundary Value Solver
// pub mod bvp;

// Tools that can be used to export data
pub mod tools {
    use polars::prelude::*;
    use std::fs::File;
    
    /// Writes a slice of column vectors to a CSV file using Polars.
    /// Each `Vec<f64>` in `data` is one column in the resulting DataFrame.
    ///
    /// # Arguments
    /// * `data` - A slice of columns, where each `Vec<i32>` is a column.
    /// * `output_path` - File path for the output CSV file.
    /// This function is to primarily export data for python plotting
    pub fn array_to_csv(
        data: &[Vec<f64>],
        cols: &[String],
        output_path: &str,
    ) -> PolarsResult<()> {
        // Basic check
        if data.is_empty() {
            return Err(PolarsError::ComputeError("Input data is empty".into()));
        }
    
        // Convert each vector into a Polars Column
        let columns: Vec<Column> = data
            .iter()
            .enumerate()
            .map(|(i, col)| {
                let col_name: PlSmallStr = cols.get(i).unwrap().clone().into(); // Convert to PlSmallStr
                Column::from(Series::new(col_name, col))
            })
            .collect();
    
        // Build the DataFrame from columns
        let mut df = DataFrame::new(columns)?;
    
        // Write out as CSV
        let file = File::create(output_path)
            .map_err(|e| PolarsError::ComputeError(format!("Could not create file: {}", e).into()))?;
    
        CsvWriter::new(file)
            .finish(&mut df)?;
    
        Ok(())
    }

    pub fn linspace(start:f64, end:f64, steps:usize) -> Vec<f64> {
        // Creates initial array and determines dx
        let mut result = Vec::with_capacity(steps);
        let dx = (end - start) / ((steps as f64) - 1.0);
    
        // Puts entries into Array
        for i in 0..steps {
            result.push(start + dx*i as f64);
        }
    
        result
    }
}



#[cfg(test)]
mod tests {
    // use super::*;



    #[test]
    fn it_works() {
        // let result = add(2, 2);
        // assert_eq!(result, 4);
    }
}
