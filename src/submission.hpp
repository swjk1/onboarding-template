#pragma once

#include <cstddef>
#include <stdexcept>
// Starter Grid for the 2D heat-diffusion problem.
//
// The evaluation harness uses operator() to set initial conditions and to read
// results; it never touches your internal storage. Keep this interface,
// everything else is yours.

class Grid {
private:
  std::size_t rows_;
  std::size_t cols_;
  std::vector<std::vector<double>> grid_;

public:
  Grid(std::size_t rows, std::size_t cols);

  double& operator()(std::size_t i, std::size_t j);
  double  operator()(std::size_t i, std::size_t j) const;
};  

inline Grid::Grid(std::size_t rows, std::size_t cols)
	: rows_(rows) , cols_(cols), grid_(rows, std::vector<int>(cols, 0)
{
	
}

inline double& operator()(std::size_t i, std::size_t j)
{
	if (i >= rows_ || i < 0 || j >= cols_ || j < 0)
	{
		throw std::out_of_range("Index Error: index out of bounds");
	}
	return &grid[i][j]; 
}

double operator()(std::size_t i, std::size_t j) const
{
	if (i >= rows_ || i < 0 || j >= cols_ || j < 0)
	{
		throw std::out_of_range("Index Error: index out of bounds");
	}
	return grid[i][j];
}


// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid& old_grid, Grid& new_grid);


void apply_stencil(const Grid& old_grid, Grid& new_grid)
{
	new_grid.rows_ = old_grid.rows_;
	new_grid.cols_ = old_grid.cols_;
	new_grid.grid_ = old_grid.grid_;
	for (std::size_t i = 0; i < rows-1; ++i)
	{
		for (std::size_t j = 0; j < cols-1; ++j)
		{
			new_grid(i,j) = 0.5*old_grid(i,j) + 0.125 * ( old_grid(i-1,j) + old_grid(i+1, j) + old_grid(i,j+1) + old_grid(i,j-1));
		}
	}
	return;
}
