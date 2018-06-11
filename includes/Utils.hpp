#include "volumetricMeshLoader.h"
#include "objMesh.h"
#include "minivector.h"
#include <string>
#include <vector>

using namespace std;

class Utils
{
	public:

		static void bouToPts(string const &vegf, vector<string> const &s_vertsf, string const &meshf, int oneIndexed=0)
		{
			VolumetricMesh * vm = VolumetricMeshLoader::load(vegf.c_str());
			ObjMesh * m = new ObjMesh(meshf);

			// m->computePseudoNormals();

			for (int i = 0; i < s_vertsf.size(); i++)
			{
				vector<int> data_vm;
				readCSV(s_vertsf[i], data_vm);
				
				vector<Vec3d> data_m;
				for (int i = 0; i < data_vm.size(); i++)
				{
					int index = data_vm[i] - 1 + oneIndexed;
					Vec3d v = vm->getVertex(index);

					int surf_mesh_index = m->getClosestVertex(v, NULL);
					Vec3d pos = m->getPosition(surf_mesh_index);
					Vec3d normal = m->getNormal(surf_mesh_index);
					Vec3d y_ref(0, 1, 0);
					Vec3d z_ref(0, 0, 1);
					Vec3d x_ref(1, 0, 0);
					
					double angle = dot(-y_ref, normal);
					if (angle >= 0 && angle <= 1)
					{
						angle = dot(-z_ref, normal);
						if (angle >= 0 && angle <= 1)
							data_m.push_back(pos);
					}
				}

				outCSV(s_vertsf[i], data_m);
			}

		}

		static void outCSV(string const &filename, vector<Vec3d> &vec)
		{
			string new_filename = filename.substr(0, filename.length() - 3) + "pts";
			ofstream outfile(new_filename.c_str());

			cout << "Writing points to: " << new_filename << endl;

			for (int i = 0; i < vec.size(); i++)
				outfile << vec[i][0] << " " << vec[i][1] << " " << vec[i][2] << endl;

			outfile.close();
		}

		static void readCSV(string const &filename, vector<int> &vec)
		{
			ifstream infile(filename.c_str());
			vector <string> record;

			while (infile)
			{
				string s;
				if (!getline( infile, s )) break;

				istringstream ss( s );

				while (ss)
				{
					string s;
					if (!getline( ss, s, ',' )) break;
						record.push_back( s );
				}
			}

			for (int i = 0; i < record.size(); i++)
			{
				vec.push_back(stoi(record[i]));
			}

			infile.close();
		}
};