#include "tutte.h"
#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
  // Load input meshes
  Eigen::MatrixXd V,U_tutte,U;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    (argc>1?argv[1]:"../shared/data/bunny.off"),V,F);
  // Load data into MatrixXd rather than VectorXd for simpler `smooth` API
  // Just use y-coordinates as data to be smoothed
  // Create a libigl Viewer object to toggle between point cloud and mesh
  igl::viewer::Viewer viewer;
  std::cout<<R"(
[space]  Toggle whether displaying 3D surface or 2D parameterization
C,c      Toggle checkerboard
t        Switch parameterization to Tutte embedding
)";
  tutte(V,F,U_tutte);

  bool plot_parameterization = false;
  const auto & update = [&]()
  {
    if(plot_parameterization)
    {
      // Viewer wants 3D coordinates, so pad UVs with column of zeros
      viewer.data.set_vertices((Eigen::MatrixXd(V.rows(),3)<<U,Eigen::VectorXd::Zero(V.rows())).finished());
    }else
    {
      viewer.data.set_vertices(V);
    }
    viewer.data.compute_normals();
    viewer.data.set_uv(U*10);
  };
  viewer.callback_key_pressed = 
    [&](igl::viewer::Viewer &, unsigned int key, int)
  {
    switch(key)
    {
      case ' ':
        plot_parameterization ^= 1;
        break;
      case 't':
        U = U_tutte;
        break;
      case 'C':
      case 'c':
        viewer.core.show_texture ^= 1;
        break;
      default:
        return false;
    }
    update();
    return true;
  };

  U = U_tutte;
  viewer.data.set_mesh(V,F);
  viewer.data.set_colors((V.array()-V.minCoeff())/(V.maxCoeff()-V.minCoeff()));
  update();
  viewer.core.show_texture = true;
  viewer.core.show_lines = false;
  viewer.launch();

  return EXIT_SUCCESS;
}
