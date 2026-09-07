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
  std::vector<std::vector<double>> grid_;

  bool index_is_valid (std::size_t i, std::size_t j);
  std::size_t get_rows();
  std::size_t get_cols();
  double& get_index(std::size_t i, std::size_t j);
  void change_rows(const std::size_t new_rows);
  void change_cols(const std::size_t new_cols);
  

  
public:
  Grid(std::size_t rows, std::size_t cols);

  double& operator()(std::size_t i, std::size_t j);
  double  operator()(std::size_t i, std::size_t j) const;
};  

inline Grid::Grid(std::size_t rows, std::size_t cols)
	: rows_(rows) , cols_(cols), grid_(rows, std::vector<double>(cols, 0.0))
{}

bool Grid::index_is_valid(std::size_t i, std::size_t j){
	if (i >= this->get_rows() || i < 0 || j >= this->get_cols() || j < 0)
	{
		return false;
	}
	return true;

}

std::size_t Grid::get_rows()
{
	return rows_;
}

std::size_t Grid::get_cols()
{
	return cols_;
}



double& Grid::get_index(std::size_t i, std::size_t j)
{
	return &grid_[i][j];
}

inline double& Grid::operator()(std::size_t i, std::size_t j)
{
	if (!index_is_valid(i,j))
	{
		throw std::out_of_range("Indicies out of range");
	}

	double& result = this->get_index(i,j);
	return result;
}

inline double Grid::operator()(std::size_t i, std::size_t j) const
{
	if (!index_is_valid(i,j))
	{
		throw std::out_of_range("Indicies out of range");
	}
	double result = get_index(i,j);
	return result;
}



// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid& old_grid, Grid& new_grid);


void apply_stencil(const Grid& old_grid, Grid& new_grid)
{
	std::size_t rows = old_grid.get_rows();
	std::size_t cols = old_grid.get_cols();
	
	for (std::size_t i = 0; i < rows; ++i)
	{
		for (std::size_t j = 0; j < cols; ++j)
		{
			if (i == 0 || j == 0 || i == rows-1 || j == cols-1)
			{
				new_grid(i,j) = old_grid(i,j);
			}
			else
			{
				new_grid(i,j) = 0.5*old_grid(i,j) + 0.125 * ( old_grid(i-1,j) + old_grid(i+1, j) + old_grid(i,j+1) + old_grid(i,j-1));
			}
		}
	}
	return;
}
