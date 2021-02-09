#include "export_fbx.h"

#include "../utils/threadutils.h"

bool export_fbx(const std::string& filename, const FBXExportData& data) {
  try {
    fbx::FBXDocument doc;
    // reverse engineered from file exported with Blender
    // doc.Objects.Geometry
    //        .Vertices
    //        .PolygonVertexIndex
    //        .LayerElementNormal (maybe unnecessary bc blender recomputes that
    //        regardless)
    //           .Normals (d)
    //        .LayerElementUV (property: i 0 or 1)
    //           .UV (d)
    //           .UVIndex (i)

    // load dummy fbx, having the correct structure for blender
    std::string dummyfile = "data/dummy.fbx";
    doc.read(dummyfile);

    fbx::FBXNode* geometry =
        doc.getNodeByName("Objects")->getChildByName("Geometry");

    // (optionally remove some broken/useless information like edges)
    geometry->removeChild("Edges");

    auto uvs = geometry->getChildrenByName("LayerElementUV");

    // naive parallel export per property, alternatively could hack fbx library
    // to set elements in parallel instead
    threadutils::parallel_for(0, 6, [&](int i) {
      if (i == 0) {
        // replace vertices
        geometry->getChildByName("Vertices")
            ->overwriteProperty(fbx::FBXProperty(data.vertices));
      } else if (i == 1) {
        // replace polygonfaceindices, where each last index of a polygon is
        // bitwise negated with ~
        geometry->getChildByName("PolygonVertexIndex")
            ->overwriteProperty(fbx::FBXProperty(data.faces));
      } else if (i == 2) {
        uvs[0]->getChildByName("UV")->overwriteProperty(
            fbx::FBXProperty(data.uv0));
      } else if (i == 3) {
        uvs[0]->getChildByName("UVIndex")->overwriteProperty(
            fbx::FBXProperty(data.uvfaces));
      } else if (i == 4) {
        uvs[1]->getChildByName("UV")->overwriteProperty(
            fbx::FBXProperty(data.uv1));
      } else {
        uvs[1]->getChildByName("UVIndex")->overwriteProperty(
            fbx::FBXProperty(data.uvfaces));
      }
    });

    // save
    doc.write(filename);
    return true;

  } catch (std::string e) {
    // std::cout << e << std::endl;
    return false;
  }

  return true;
}