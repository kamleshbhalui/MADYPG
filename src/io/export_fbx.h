#include "../dependencies/fbx/fbxdocument.h"

struct FBXExportData {
  std::vector<float> vertices; // flat x0 y0 z0 x1 y1 z1 ...
  std::vector<int32_t> faces; // flat v0 v1 ~v2 v3 v4 ~v5 ... last poly index negated with ~
  std::vector<float> uv0; // flat u0 v0 u1 v1 ... first uvmap
  std::vector<float> uv1; // flat u0 v0 u1 v1 ... second uvmap
  std::vector<int32_t> uvfaces; // flat v0 v1 v2 v3 v4 v5 ... (unnegated!) assumed same for both uv maps
};

bool export_fbx(const std::string& filename, const FBXExportData& data);