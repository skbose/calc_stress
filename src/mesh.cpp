#include "../includes/mesh.hpp"

void triangle_mesh_t::clear()
{
  for (auto it = triangles.begin(); it != triangles.end(); ++it)
  {
    delete *it;
  }

  for (auto it = vertices.begin(); it != vertices.end(); ++it)
  {
    delete *it;
  }
}

bool get_first_value(std::istringstream & ss, int & value)
{
  std::string s;

  if (ss >> s)
  {
    std::istringstream valuestream(s);
    if (valuestream >> value)
      return true;
  }

  return false;
}

bool get_first_value(std::istringstream & ss, double & value)
{
  std::string s;

  if (ss >> s)
  {
    std::istringstream valuestream(s);
    if (valuestream >> value)
      return true;
  }

  return false;
}

bool triangle_mesh_t::load(std::string const & path)
{
  clear();

  std::ifstream f(path.c_str());
  if (!f)
  {
    std::cerr << "Mesh: Couldn't load file " << path << std::endl;
    return false;
  }

  std::string line;
  // int normals = 0;

  while (std::getline(f, line))
  {
    std::istringstream linestream(line);
    std::string op;
    linestream >> op;

    if (op[0] == '#')
      continue;

    if (op == "v")
    {
      double x, y, z;

      linestream >> x >> y >> z;
      Eigen::Vector3d v(x, y, z);
      
      vertices.push_back(new mesh_vertex(v));
    }

    if (op == "f")   // extract a face as a triangle fan
    {
      int first, second, third;

      if (!get_first_value(linestream, first))
        continue;

      if (!get_first_value(linestream, second))
        continue;

      while (get_first_value(linestream, third))
      {
        Eigen::Vector3d p1 = vertices[first - 1]->pos;
        Eigen::Vector3d p2 = vertices[second - 1]->pos;
        Eigen::Vector3d p3 = vertices[third - 1]->pos;

        Eigen::Vector3d p12 = p2 - p1;
        Eigen::Vector3d p13 = p3 - p1;

        double area_ = (p12.cross(p13)).norm() / 2.0;

        triangles.push_back(new mesh_triangle(first - 1, second - 1, third - 1, area_));
        second = third;
      }
    }
  }

  return true;
}

void triangle_mesh_t::flushVertexData()
{
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i]->area = 0.0;
    vertices[i]->ratio = 0.0;
    vertices[i]->neighborAreaDeformation = 0.0;
    vertices[i]->neighborAreaOriginal = 0.0;
    vertices[i]->lx = vertices[i]->ly = vertices[i]->lz = 0.0;
  
    vertices[i]->neighbours = 0;
  }  
}

void triangle_mesh_t::calculatePtArea()
{
  for (auto it = triangles.begin(); it != triangles.end(); ++it)
  {
    int v_index_0 = (*it)->v_ind[0];
    int v_index_1 = (*it)->v_ind[1];
    int v_index_2 = (*it)->v_ind[2];

    vertices[v_index_0]->area += 0.33 * (*it)->area;
    vertices[v_index_1]->area += 0.33 * (*it)->area;
    vertices[v_index_2]->area += 0.33 * (*it)->area;
  }
}

bool triangle_mesh_t::savePtArea(std::string const & path)
{
  std::ofstream f(path.c_str());
  
  if (!f)
  {
    std::cerr << "Couldn't open file " << path << std::endl;
    return false;
  }

  for (auto it = vertices.begin(); it != vertices.end(); ++it)
  {
    f << (*it)->area << std::endl;  
  }

  f.close();

  return true;
}

bool triangle_mesh_t::saveMassRatio(std::string const & in_file_path, std::string const & out_file_path)
{
  std::fstream f(in_file_path.c_str(), std::fstream::in);
  std::vector<double> sdf_values;

  if (!f)
  {
    std::cerr << "Couldn't open file " << in_file_path << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(f, line))
  {
    std::istringstream linestream(line);
    double v1, v2, v3, sdf;

    linestream >> v1 >> v2 >> v3 >> sdf;
    sdf_values.push_back(sdf);
  }

  double denom = 0.0;
  for (int i = 0; i < vertices.size(); i++)
  {
    denom += vertices[i]->area * sdf_values[i];
  }

  f.close();

  f.open(out_file_path.c_str(), std::fstream::out);
  
  if (!f)
  {
    std::cerr << "Couldn't open file " << out_file_path << std::endl;
    return false;
  }

  for (int i = 0; i < vertices.size(); i++)
  {
    double ratio_ = (vertices[i]->area * sdf_values[i]) / denom;
    vertices[i]->ratio = ratio_;
    
    f << ratio_ << std::endl;
  }

  f.close();

  return true;
}

bool triangle_mesh_t::saveMassRatioAsFeature(std::string const & path)
{
  std::ofstream f(path.c_str());
  
  if (!f)
  {
    std::cerr << "Couldn't open file " << path << std::endl;
    return false;
  }

  for (auto it = vertices.begin(); it != vertices.end(); ++it)
  {
    double x = (*it)->pos.x();
    double y = (*it)->pos.y();
    double z = (*it)->pos.z();

    f << x << " " << y << " " << z << " " << (*it)->ratio << std::endl;  
  }

  f.close();

  return true;
}

bool triangle_mesh_t::loadDeformation(std::string const & path)
{
  std::ifstream f(path);
  
  if (!f)
  {
    std::cerr << "Couldn't open file " << path << std::endl;
    return false;
  }

  flushVertexData();

  int index = 0;
  std::string line;

  while (getline(f, line))
  {
    std::istringstream ss(line);
    double ux, uy, uz;

    ss >> ux >> uy >> uz;

    vertices[index]->ux = ux;
    vertices[index]->uy = uy;
    vertices[index]->uz = uz;

    index++;
  }

  f.close();

  return true;
}

bool triangle_mesh_t::loadDeformation(double * u)
{
  flushVertexData();
  
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i]->ux = u[3*i+0];
    vertices[i]->uy = u[3*i+1];
    vertices[i]->uz = u[3*i+2];    
  }

  return true;
}

bool triangle_mesh_t::laplacianOnDeformations()
{
  for (auto it = triangles.begin(); it != triangles.end(); ++it)
  {
    int v_index_0 = (*it)->v_ind[0];
    int v_index_1 = (*it)->v_ind[1];
    int v_index_2 = (*it)->v_ind[2];

    vertices[v_index_0]->lx += vertices[v_index_1]->ux + vertices[v_index_2]->ux;
    vertices[v_index_0]->ly += vertices[v_index_1]->uy + vertices[v_index_2]->uy;
    vertices[v_index_0]->lz += vertices[v_index_1]->uz + vertices[v_index_2]->uz;
    vertices[v_index_0]->neighbours += 2;

    vertices[v_index_1]->lx += vertices[v_index_2]->ux + vertices[v_index_0]->ux;
    vertices[v_index_1]->ly += vertices[v_index_2]->uy + vertices[v_index_0]->uy;
    vertices[v_index_1]->lz += vertices[v_index_2]->uz + vertices[v_index_0]->uz;
    vertices[v_index_1]->neighbours += 2;

    vertices[v_index_2]->lx += vertices[v_index_1]->ux + vertices[v_index_0]->ux;
    vertices[v_index_2]->ly += vertices[v_index_1]->uy + vertices[v_index_0]->uy;
    vertices[v_index_2]->lz += vertices[v_index_1]->uz + vertices[v_index_0]->uz;
    vertices[v_index_2]->neighbours += 2;
  }

  for (auto it = vertices.begin(); it != vertices.end(); ++it)
  {
    (*it)->lx /= 2.0;
    (*it)->neighbours /= 2;
    (*it)->lx -= (*it)->neighbours * (*it)->ux;

    (*it)->ly /= 2.0;
    (*it)->neighbours /= 2;
    (*it)->ly -= (*it)->neighbours * (*it)->uy;

    (*it)->lz /= 2.0;
    (*it)->neighbours /= 2;
    (*it)->lz -= (*it)->neighbours * (*it)->uz;
  }

  return true;
}

bool triangle_mesh_t::saveAbsLaplacianValuesAsFeature(std::string const & path)
{

  std::ofstream f(path);
  
  if (!f)
  {
    std::cerr << "Couldn't open file " << path << std::endl;
    return false;
  }

  for (auto it = vertices.begin(); it != vertices.end(); ++it)
  {
    double x = (*it)->pos.x();
    double y = (*it)->pos.y();
    double z = (*it)->pos.z();

    double lx = fabs((*it)->lx);
    double ly = fabs((*it)->ly);
    double lz = fabs((*it)->lz);

    f << x << " " << y << " " << z << " " << lx << " " << ly << " " << lz << std::endl; 
  }

  f.close();
}

bool triangle_mesh_t::saveLaplacianVectorNormAsFeature(std::string const & path, int max_points)
{
  std::ofstream f(path);
  
  if (!f)
  {
    std::cerr << "Couldn't open file " << path << std::endl;
    return false;
  }

  int num_points = 0;
  bool disp_all = max_points == -1 ? true : false;
  for (auto it = vertices.begin(); (num_points < max_points || disp_all) && it != vertices.end(); ++it, num_points++)
  {
    double x = (*it)->pos.x();
    double y = (*it)->pos.y();
    double z = (*it)->pos.z();

    double lx = fabs((*it)->lx);
    double ly = fabs((*it)->ly);
    double lz = fabs((*it)->lz);

    Eigen::Vector3d l(lx, ly, lz);

    f << x << " " << y << " " << z << " " << l.norm() << std::endl; 
  }

  f.close(); 
}

void triangle_mesh_t::calculateLocalAreaDeformationAboutPoints()
{
  // do the pre-processing first
	calculateAreaDeformationOfMeshTriangles();

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
  {
    int v_index_0 = (*it)->v_ind[0];
    int v_index_1 = (*it)->v_ind[1];
    int v_index_2 = (*it)->v_ind[2];

    // double abs_area_deformation = fabs((*it)->area - (*it)->deformedArea);
    double area_original =  (*it)->area;
    double area_deformed = (*it)->deformedArea;

    double ratio1 = area_original / area_deformed;
    double ratio2 = area_deformed / area_original;
    double distort_metric = fmax(fabs(ratio1), fabs(ratio2)) - 1;

    vertices[v_index_0]->neighborAreaDeformation += fabs(distort_metric) * area_original;
    vertices[v_index_1]->neighborAreaDeformation += fabs(distort_metric) * area_original;
    vertices[v_index_2]->neighborAreaDeformation += fabs(distort_metric) * area_original;

    vertices[v_index_0]->neighborAreaOriginal += area_original;
    vertices[v_index_1]->neighborAreaOriginal += area_original;
    vertices[v_index_2]->neighborAreaOriginal += area_original;
  }

  for (auto it = vertices.begin(); it != vertices.end(); ++it)
  {
  	// avg area deformation captured here.
    if ((*it)->neighborAreaOriginal > 1e-10)
      (*it)->neighborAreaDeformation /= (*it)->neighborAreaOriginal;
  }

}

void triangle_mesh_t::calculateAreaDeformationOfMeshTriangles()
{
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
	  int v_index_0 = (*it)->v_ind[0];
  	int v_index_1 = (*it)->v_ind[1];
  	int v_index_2 = (*it)->v_ind[2];

    // calculate deformation vector of the first vertex
    double ux_ = vertices[v_index_0]->ux;
    double uy_ = vertices[v_index_0]->uy;
    double uz_ = vertices[v_index_0]->uz;
    Eigen::Vector3d d_p0(ux_, uy_, uz_);
    
    // calculate deformation vector of the second vertex
    ux_ = vertices[v_index_1]->ux;
    uy_ = vertices[v_index_1]->uy;
    uz_ = vertices[v_index_1]->uz;
    Eigen::Vector3d d_p1(ux_, uy_, uz_);

    // calculate deformation vector of the third vertex
    ux_ = vertices[v_index_2]->ux;
    uy_ = vertices[v_index_2]->uy;
    uz_ = vertices[v_index_2]->uz;
    Eigen::Vector3d d_p2(ux_, uy_, uz_);


    Eigen::Vector3d p1 = vertices[v_index_0]->pos + d_p0;
    Eigen::Vector3d p2 = vertices[v_index_1]->pos + d_p1;
    Eigen::Vector3d p3 = vertices[v_index_2]->pos + d_p2;

    Eigen::Vector3d p12 = p2 - p1;
    Eigen::Vector3d p13 = p3 - p1;

    double area_new = (p12.cross(p13)).norm() / 2.0;

    (*it)->deformedArea = area_new;
	}
}

bool triangle_mesh_t::saveLocalAreaDeformationAboutPointsAsFeature(std::string const &path, int max_points)
{
  std::ofstream f(path);
  
  if (!f)
  {
    std::cerr << "Couldn't open file " << path << std::endl;
    return false;
  }

  int num_points = 0;
  bool disp_all = max_points == -1 ? true : false;
  for (auto it = vertices.begin(); (num_points < max_points || disp_all) && it != vertices.end(); ++it, num_points++)
  {
    double x = (*it)->pos.x();
    double y = (*it)->pos.y();
    double z = (*it)->pos.z();

    f << x << " " << y << " " << z << " " << (*it)->neighborAreaDeformation << std::endl; 
  }

  f.close();

  std::cout << "Area distortion written to: " << path << std::endl; 
}

bool triangle_mesh_t::savePointMovementMagnitudeAsFeature(std::string const &path, int max_points)
{
	std::ofstream f(path);
  
	if (!f)
	{
	std::cerr << "Couldn't open file " << path << std::endl;
	return false;
	}

	int num_points = 0;
	bool disp_all = max_points == -1 ? true : false;
	for (int i = 0; (num_points < max_points || disp_all) && i < this->vertices.size(); i++, num_points++)
	{
		double x = (this)->vertices[i]->pos.x();
    double y = (this)->vertices[i]->pos.y();
    double z = (this)->vertices[i]->pos.z();

    double ux = (this->vertices[i])->ux;
    double uy = (this->vertices[i])->uy;
    double uz = (this->vertices[i])->uz;

		Eigen::Vector3d movementDir(ux, uy, uz);
		double movementDirNorm = movementDir.norm();

		f << x << " " << y << " " << z << " " << movementDirNorm << std::endl; 
	}

	f.close();

  std::cout << "Points displacement written to: " << path << std::endl;
}

void triangle_mesh_t::operator-(triangle_mesh_t const &m)
{
  for (int i = 0; i < this->vertices.size(); i++)
  {
    double x = (this)->vertices[i]->pos.x();
    double y = (this)->vertices[i]->pos.y();
    double z = (this)->vertices[i]->pos.z();

    double ux = this->vertices[i]->pos.x() - m.vertices[i]->pos.x();
    double uy = this->vertices[i]->pos.y() - m.vertices[i]->pos.y();
    double uz = this->vertices[i]->pos.z() - m.vertices[i]->pos.z();

    Eigen::Vector3d movementDir(ux, uy, uz);
    double movementDirNorm = movementDir.norm();
    
    std::cout << -ux << " " << -uy << " " << -uz << std::endl; 
  }
}

bool triangle_mesh_t::savePointAndDeformation(std::string const &path, int max_points)
{
  std::ofstream f(path);
  
  if (!f)
  {
  std::cerr << "Couldn't open file " << path << std::endl;
  return false;
  }

  int num_points = 0;
  bool disp_all = max_points == -1 ? true : false;
  for (int i = 0; (num_points < max_points || disp_all) && i < this->vertices.size(); i++, num_points++)
  {
    double x = (this)->vertices[i]->pos.x();
    double y = (this)->vertices[i]->pos.y();
    double z = (this)->vertices[i]->pos.z();

    double ux = (this->vertices[i])->ux;
    double uy = (this->vertices[i])->uy;
    double uz = (this->vertices[i])->uz;

    f << x << " " << y << " " << z << " " << ux << " " << uz << " " << uy << std::endl; 
  }

  f.close(); 
  
  return true;
}