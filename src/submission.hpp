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

  bool index_is_valid (std::size_t i, std::size_t j) const;
  std::size_t get_rows() const;
  std::size_t get_cols() const;
  double& get_index(std::size_t i, std::size_t j);
  double const_get_index(std::size_t i, std::size_t j) const;
  void change_rows(const std::size_t new_rows);
  void change_cols(const std::size_t new_cols);
  
  

  
public:
  Grid(std::size_t rows, std::size_t cols);

  double& operator()(std::size_t i, std::size_t j);
  double  operator()(std::size_t i, std::size_t j) const;
  const std::size_t* get_dims() const;
};  

inline const std::size_t* Grid::get_dims() const
{
	std::size_t rows = this->get_rows();
	std::size_t cols = this->get_cols();
	const std::size_t* dims = new std::size_t [2] {rows, cols};
	return dims;
}

inline Grid::Grid(std::size_t rows, std::size_t cols)
	: rows_(rows) , cols_(cols), grid_(rows * cols)
{}

inline bool Grid::index_is_valid(std::size_t i, std::size_t j) const{
	if (i >= this->get_rows() || j >= this->get_cols())
	{
		return false;
	}
	return true;

}

inline std::size_t Grid::get_rows() const
{
	return rows_;
}

inline std::size_t Grid::get_cols() const
{
	return cols_;
}



inline double& Grid::get_index(std::size_t i, std::size_t j)
{
	return grid_[i*(this->get_cols()) + j];
}

inline double Grid::const_get_index(std::size_t i, std::size_t j) const
{
	return grid_[(i*this->get_cols()) + j];
}

inline double& Grid::operator()(std::size_t i, std::size_t j)
{

	double& result = this->get_index(i,j);
	return result;
}

inline double Grid::operator()(std::size_t i, std::size_t j) const
{
	const double result = this->const_get_index(i,j);
	return result;
}



// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid& old_grid, Grid& new_grid);


void apply_stencil(const Grid& old_grid, Grid& new_grid)
{	
	const std::size_t* dims = old_grid.get_dims();
	std::size_t rows = *(dims);
	std::size_t cols = *(dims + 1);
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
	delete[] dims;
	dims=nullptr;
	return;
}
