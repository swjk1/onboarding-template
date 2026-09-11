#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>
#include <array>
// Starter Grid for the 2D heat-diffusion problem.
//
// The evaluation harness uses operator() to set initial conditions and to read
// results; it never touches your internal storage. Keep this interface,
// everything else is yours.

class Grid {
private:
  std::size_t rows_;
  std::size_t cols_;
  std::vector<double> grid_;

  double& get_index(std::size_t i, std::size_t j);
  const double& get_index(std::size_t i, std::size_t j) const;
  

  
public:
  Grid(std::size_t rows, std::size_t cols);
  
  double& operator()(std::size_t i, const std::size_t j);
  double  operator()(std::size_t i, std::size_t j) const;
  std::array<std::size_t, 2> get_dims() const;
  const double* get_data() const;
  double* get_data();
  
};



const bool validate(const Grid& grid1, const Grid& grid2);

const double* Grid::get_data() const
{
	return grid_.data();
}

double* Grid::get_data()
{
	return grid_.data();
}


inline std::array<std::size_t, 2> Grid::get_dims() const
{
	return {rows_, cols_};
}

inline Grid::Grid(std::size_t rows, std::size_t cols)
	: rows_(rows) , cols_(cols), grid_(rows * cols)
{}



inline double& Grid::get_index(std::size_t i, std::size_t j)
{
	return grid_[i*cols_ + j];
}

inline const double& Grid::get_index(std::size_t i, std::size_t j) const
{
	return grid_[i*cols_ + j];
}

inline double& Grid::operator()(std::size_t i, std::size_t j)
{

	return get_index(i,j);
}

inline double Grid::operator()(std::size_t i, std::size_t j) const
{
	return get_index(i,j);
}


inline const bool validate(const Grid& grid1, const Grid& grid2)
{
	return grid1.get_dims() == grid2.get_dims();
}
// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid& old_grid, Grid& new_grid);


void apply_stencil(const Grid& old_grid, Grid& new_grid)
{	
	const std::array<std::size_t, 2>& dims = old_grid.get_dims();
	const std::size_t& rows = dims[0];
	const std::size_t& cols = dims[1];
	const double* old_grid_data = old_grid.get_data();
	double* new_grid_data = new_grid.get_data();
	for (std::size_t i {1}; i < rows-1; ++i)
	{
		for (std::size_t k {i*cols + 1}; k < i*cols + cols-1; ++k)
		{
			new_grid_data[k] = 0.5*old_grid_data[k] + 0.125 * ( old_grid_data[k-cols] + old_grid_data[k+cols] + old_grid_data[k+1] + old_grid_data[k-1]);
		}
	}
	std::size_t c1 { cols-1 };
	std::size_t r1 { rows-1 };
	for (std::size_t i {}; i < rows; ++i)
	{
		new_grid(i,0) = old_grid(i,0);
		new_grid(i,c1) = old_grid(i,c1);
	}
	for (std::size_t j {}; j < cols; ++j)
	{
		new_grid(0,j) = old_grid(0,j);
		new_grid(r1,j) = old_grid(r1, j);
	}
	return;
}
