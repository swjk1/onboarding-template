#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>
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

  std::size_t get_rows() const;
  std::size_t get_cols() const;
  double& get_index(const std::size_t& i, const std::size_t& j);
  const double& get_index(const std::size_t& i, const std::size_t& j) const;
  
  

  
public:
  Grid(std::size_t rows, std::size_t cols);
  
  double& operator()(const std::size_t& i, const std::size_t& j);
  double  operator()(const std::size_t& i, const std::size_t& j) const;
  std::vector<std::size_t> get_dims() const;
  
  
};

const bool validate(const Grid& grid1, const Grid& grid2);
  

inline std::vector<std::size_t> Grid::get_dims() const
{
	const std::vector<std::size_t> dims {get_rows(), get_cols()};
	return dims;
}

inline Grid::Grid(std::size_t rows, std::size_t cols)
	: rows_(rows) , cols_(cols), grid_(rows * cols)
{}

inline std::size_t Grid::get_rows() const
{
	return rows_;
}

inline std::size_t Grid::get_cols() const
{
	return cols_;
}



inline double& Grid::get_index(const std::size_t& i, const std::size_t& j)
{
	return grid_[i*(this->get_cols()) + j];
}

inline const double& Grid::get_index(const std::size_t& i, const std::size_t& j) const
{
	return grid_[(i*this->get_cols()) + j];
}

inline double& Grid::operator()(const std::size_t& i, const std::size_t& j)
{

	return get_index(i,j);
}

inline double Grid::operator()(const std::size_t& i, const std::size_t& j) const
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
	const std::vector<std::size_t>& dims = old_grid.get_dims();
	const std::size_t& rows = dims[0];
	const std::size_t& cols = dims[1];
	for (std::size_t i {1}; i < rows; ++i)
	{
		for (std::size_t j {1}; j < cols; ++j)
		{
			new_grid(i,j) = 0.5*old_grid(i,j) + 0.125 * ( old_grid(i-1,j) + old_grid(i+1, j) + old_grid(i,j+1) + old_grid(i,j-1));
		}
	}
	for (std::size_t i {}; i < rows; ++i)
	{
		new_grid(i,0) = old_grid(i,0);
		new_grid(i,cols-1) = old_grid(i,cols-1);
	}
	for (std::size_t j {}; j < cols; ++j)
	{
		new_grid(0,j) = old_grid(0,j);
		new_grid(rows-1,j) = old_grid(rows-1, j);
	}
	return;
}
