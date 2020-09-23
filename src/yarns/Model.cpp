
#include "Model.h"


#include <fstream>
#include "../utils/debug_logging.h"

void deserialize_matrix(const std::string& filename, MatrixGLf &M) {
  std::ifstream ifs(filename);
  if (!ifs) {
    Debug::msgassert("Couldn't load matrix file!", bool(ifs));
    return;
  }

  std::string line;
  int i = 0;
  while (std::getline(ifs, line)) {
    std::stringstream ss(line);
    if(i >= M.rows()) {
      Debug::msgassert("Matrix needs to be allocated before loading!", false);
      return;
    }
    for (int j = 0; j < M.cols(); j++) {
      ss >> M(i, j);
    }
    ++i;
  }
}

const std::vector<std::string> Model::strain_names = {"sx","sa","sy","IIx","IIy"};