#define _USE_MATH_DEFINES
#include <cmath>
#include "iostream"
#include <string>
#include "patchOrganizerS.h"
#include "findMatch.h"
#include "bezier.h"

using namespace PMVS3;
using namespace Patch;
using namespace std;

Ppatch CpatchOrganizerS::m_MAXDEPTH(new Cpatch());
Ppatch CpatchOrganizerS::m_BACKGROUND(new Cpatch());

CpatchOrganizerS::CpatchOrganizerS(CfindMatch& findMatch) : m_fm(findMatch) {
}

// change the contents of m_images from images to indexes
void CpatchOrganizerS::image2index(Cpatch& patch) {
  // first image has to be target imag
  vector<int> newimages;
  for (int i = 0; i < (int)patch.m_images.size(); ++i) {
    const int index = m_fm.m_pss.image2index(patch.m_images[i]);
    if (index != -1)
      newimages.push_back(index);
  }
  
  patch.m_images.swap(newimages);  

  // make sure that the reference image is the targeting image
  int exist = -1;
  for (int j = 0; j < (int)patch.m_images.size(); ++j) {
    if (patch.m_images[j] < m_fm.m_tnum) {
      exist = j;
      break;
    }
  }
  if (exist == -1)
    patch.m_images.clear();
  else if (exist != 0)
    swap(patch.m_images[0], patch.m_images[exist]);
}

// change the contents of m_images from indexes to images
void CpatchOrganizerS::index2image(Cpatch& patch) {
  for (int i = 0; i < (int)patch.m_images.size(); ++i)
    patch.m_images[i] = m_fm.m_pss.m_images[patch.m_images[i]];
  for (int i = 0; i < (int)patch.m_vimages.size(); ++i)
    patch.m_vimages[i] = m_fm.m_pss.m_images[patch.m_vimages[i]];
}

void CpatchOrganizerS::init(void) {
  m_pgrids.clear();   m_pgrids.resize(m_fm.m_tnum);
  m_vpgrids.clear();  m_vpgrids.resize(m_fm.m_tnum);
  m_dpgrids.clear();  m_dpgrids.resize(m_fm.m_tnum);
  m_counts.clear();   m_counts.resize(m_fm.m_tnum);
  
  m_gwidths.clear();  m_gwidths.resize(m_fm.m_num);
  m_gheights.clear(); m_gheights.resize(m_fm.m_num);
  for (int index = 0; index < m_fm.m_num; ++index) {
    const int gwidth = (m_fm.m_pss.getWidth(index, m_fm.m_level)
                        + m_fm.m_csize - 1) / m_fm.m_csize;
    const int gheight = (m_fm.m_pss.getHeight(index, m_fm.m_level)
                         + m_fm.m_csize - 1) / m_fm.m_csize;
    m_gwidths[index] = gwidth;
    m_gheights[index] = gheight;

    if (index < m_fm.m_tnum) {
      m_pgrids[index].resize(gwidth * gheight);
      m_vpgrids[index].resize(gwidth * gheight);
      m_dpgrids[index].resize(gwidth * gheight);
      m_counts[index].resize(gwidth * gheight);
      fill(m_dpgrids[index].begin(), m_dpgrids[index].end(), m_MAXDEPTH);
    }
  }
}

void CpatchOrganizerS::writePatches2(const std::string prefix, bool bExportPLY, bool bExportPatch, bool bExportPSet) {
  collectPatches(1);
  
  if (bExportPLY)
  {
    char buffer[1024];
    sprintf(buffer, "%s.ply", prefix.c_str());
    writePLY(m_ppatches, buffer);
  }

  if (bExportPatch)
  {
    char buffer[1024];
    sprintf(buffer, "%s.patch", prefix.c_str());
    ofstream ofstr;
    ofstr.open(buffer);
    ofstr << "PATCHES" << endl
          << (int)m_ppatches.size() << endl;
    for (int p = 0; p < (int)m_ppatches.size(); ++p) {
      Cpatch patch = *m_ppatches[p];
      index2image(patch);
      ofstr << patch << "\n";
    }
    ofstr.close();
  }

  if (bExportPSet)
  {
    char buffer[1024];
    sprintf(buffer, "%s.pset", prefix.c_str());
    ofstream ofstr;
    ofstr.open(buffer);
    for (int p = 0; p < (int)m_ppatches.size(); ++p)
      ofstr << m_ppatches[p]->m_coord[0] << ' '
            << m_ppatches[p]->m_coord[1] << ' '
            << m_ppatches[p]->m_coord[2] << ' '
            << m_ppatches[p]->m_normal[0] << ' '
            << m_ppatches[p]->m_normal[1] << ' '
            << m_ppatches[p]->m_normal[2] << "\n";
    ofstr.close();
  }
}

void CpatchOrganizerS::readPatches(void) {
  // Read-in existing reconstructed points. set m_fix to one for non-targeting images
  for (int i = 0; i < m_fm.m_tnum; ++i) {
    const int image = m_fm.m_images[i];
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.patc%d",
            m_fm.m_prefix.c_str(), image, m_fm.m_level);
    ifstream ifstr;
    ifstr.open(buffer);
    if (!ifstr.is_open())
      continue;

    string header;    int pnum;
    ifstr >> header >> pnum;
    cerr << image << ' ' << pnum << " patches" << endl;
    for (int p = 0; p < pnum; ++p) {
      Ppatch ppatch(new Cpatch());
      ifstr >> *ppatch;
      ppatch->m_fix = 0;
      ppatch->m_vimages.clear();

      image2index(*ppatch);
      if (ppatch->m_images.empty())
        continue;
      
      // m_vimages must be targeting images
#ifdef DEBUG
      for (int j = 0; j < (int)ppatch->m_vimages.size(); ++j)
        if (m_fm.m_tnum <= ppatch->m_vimages[j]) {
          cerr << "Impossible in readPatches. m_vimages must be targeting images" << endl
               << "for patches stored in targeting images, if visdata2 have been consistent" << endl;
          exit (1);
        }
#endif
      setGrids(*ppatch);
      addPatch(ppatch);
    }
    
    ifstr.close();
  }

  // For patches in non-targeting images
  for (int i = m_fm.m_tnum; i < m_fm.m_num; ++i) {
    const int image = m_fm.m_images[i];
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.patc%d", m_fm.m_prefix.c_str(), image, m_fm.m_level);
    ifstream ifstr;
    ifstr.open(buffer);
    if (!ifstr.is_open())
      continue;
    
    string header;    int pnum;
    ifstr >> header >> pnum;
    cerr << image << ' ' << pnum << " patches" << endl;
    
    for (int p = 0; p < pnum; ++p) {
      Ppatch ppatch(new Cpatch());
      ifstr >> *ppatch;      
      ppatch->m_fix = 1;
      ppatch->m_vimages.clear();
      
      image2index(*ppatch);
      if (ppatch->m_images.empty())
        continue;
      
      setGrids(*ppatch);
      addPatch(ppatch);
    }
    
    ifstr.close();
  }
}

void CpatchOrganizerS::collectPatches(const int target) {
  m_ppatches.clear();

  for (int index = 0; index < m_fm.m_tnum; ++index) {
    for (int i = 0; i < (int)m_pgrids[index].size(); ++i) {
      vector<Ppatch>::iterator begin = m_pgrids[index][i].begin();
      while (begin != m_pgrids[index][i].end()) {
        (*begin)->m_id = -1;
        begin++;
      }
    }
  }
  
  int count = 0;
  for (int index = 0; index < m_fm.m_tnum; ++index) {
    for (int i = 0; i < (int)m_pgrids[index].size(); ++i) {
      vector<Ppatch>::iterator begin = m_pgrids[index][i].begin();
      while (begin != m_pgrids[index][i].end()) {
        if ((*begin)->m_id == -1) {
          (*begin)->m_id = count++;

          if (target == 0 || (*begin)->m_fix == 0)
            m_ppatches.push_back(*begin);
        }
        ++begin;
      }
    }
  }
}

void CpatchOrganizerS::collectPatches(std::priority_queue<Patch::Ppatch,
                                      std::vector<Patch::Ppatch>,
                                      P_compare>& pqpatches) {
  for (int index = 0; index < m_fm.m_tnum; ++index) {
    for (int i = 0; i < (int)m_pgrids[index].size(); ++i) {
      vector<Ppatch>::iterator begin = m_pgrids[index][i].begin();
      while (begin != m_pgrids[index][i].end()) {
        if ((*begin)->m_flag == 0) {
          (*begin)->m_flag = 1;
          pqpatches.push(*begin);
        }
        ++begin;
      }
    }
  }
}

void CpatchOrganizerS::collectPatches(const int index,
                                     std::priority_queue<Patch::Ppatch, std::vector<Patch::Ppatch>,
                                      P_compare>& pqpatches) {
  m_fm.m_imageLocks[index].wrlock();
  for (int i = 0; i < (int)m_pgrids[index].size(); ++i) {
    vector<Ppatch>::iterator begin = m_pgrids[index][i].begin();
    vector<Ppatch>::iterator end = m_pgrids[index][i].end();
    
    while (begin != end) {
      if ((*begin)->m_images[0] == index && (*begin)->m_flag == 0) {
        (*begin)->m_flag = 1;
        pqpatches.push(*begin);
      }
      ++begin;
    }
  }
  m_fm.m_imageLocks[index].unlock();
}

// Should be used only for writing
void CpatchOrganizerS::collectNonFixPatches(const int index,
                                      std::vector<Patch::Ppatch>& ppatches) {
  m_fm.m_imageLocks[index].wrlock();;
  for (int i = 0; i < (int)m_pgrids[index].size(); ++i) {
    vector<Ppatch>::iterator begin = m_pgrids[index][i].begin();
    vector<Ppatch>::iterator end = m_pgrids[index][i].end();
    
    while (begin != end) {
      if ((*begin)->m_images[0] == index && (*begin)->m_fix == 0) {
        ppatches.push_back(*begin);
      }
      ++begin;
    }
  }
  m_fm.m_imageLocks[index].unlock();
}

void CpatchOrganizerS::clearFlags(void) {
  vector<Ppatch>::iterator bppatch = m_ppatches.begin();
  vector<Ppatch>::iterator eppatch = m_ppatches.end();

  while (bppatch != eppatch) {
    (*bppatch)->m_flag = 0;
    ++bppatch;
  }
}

void CpatchOrganizerS::clearCounts(void) {
  for (int index = 0; index < m_fm.m_tnum; ++index) {
    vector<unsigned char>::iterator begin = m_counts[index].begin();
    vector<unsigned char>::iterator end = m_counts[index].end();
    while (begin != end) {
      *begin = (unsigned char)0;
      ++begin;
    }
  }
}

void CpatchOrganizerS::addPatch(Patch::Ppatch& ppatch) {
  // First handle m_vimages
  vector<int>::iterator bimage = ppatch->m_images.begin();
  vector<int>::iterator eimage = ppatch->m_images.end();
  vector<Vec2i>::iterator bgrid = ppatch->m_grids.begin();
  while (bimage != eimage) {
    const int index = *bimage;
    if (m_fm.m_tnum <= index) {
      ++bimage;      ++bgrid;
      continue;
    }
    
    const int index2 = (*bgrid)[1] * m_gwidths[index] + (*bgrid)[0];
    m_fm.m_imageLocks[index].wrlock();
    m_pgrids[index][index2].push_back(ppatch);
	m_fm.m_imageLocks[index].unlock();
    ++bimage;
    ++bgrid;
  }

  // If depth, set vimages
  if (m_fm.m_depth == 0)
    return;
      
  bimage = ppatch->m_vimages.begin();
  eimage = ppatch->m_vimages.end();
  bgrid = ppatch->m_vgrids.begin();
  
  while (bimage != eimage) {
    const int index = *bimage;
    const int index2 = (*bgrid)[1] * m_gwidths[index] + (*bgrid)[0];
	m_fm.m_imageLocks[index].wrlock();
    m_vpgrids[index][index2].push_back(ppatch);
	m_fm.m_imageLocks[index].unlock();
    ++bimage;
    ++bgrid;
  }

  updateDepthMaps(ppatch);
}

void CpatchOrganizerS::updateDepthMaps(Ppatch& ppatch) {
  for (int image = 0; image < m_fm.m_tnum; ++image) {
    const Vec3f icoord = m_fm.m_pss.project(image, ppatch->m_coord, m_fm.m_level);

    const float fx = icoord[0] / m_fm.m_csize;
    const int xs[2] = {(int)floor(fx), (int)ceil(fx)};
    const float fy = icoord[1] / m_fm.m_csize;
    const int ys[2] = {(int)floor(fy), (int)ceil(fy)};
    
    const float depth = m_fm.m_pss.m_photos[image].m_oaxis * ppatch->m_coord;

	m_fm.m_imageLocks[image].wrlock();
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
	if (xs[i] < 0 || m_gwidths[image] <= xs[i] ||
            ys[j] < 0 || m_gheights[image] <= ys[j])
	  continue;

        const int index = ys[j] * m_gwidths[image] + xs[i];
	if (m_dpgrids[image][index] == m_MAXDEPTH)
	  m_dpgrids[image][index] = ppatch;
	else {
	  const float dtmp = m_fm.m_pss.m_photos[image].m_oaxis *
	    m_dpgrids[image][index]->m_coord;
	  
	  if (depth < dtmp)
	    m_dpgrids[image][index] = ppatch;
	}
      }
    }
    m_fm.m_imageLocks[image].unlock();
  }
}

void CpatchOrganizerS::setGridsImages(Patch::Cpatch& patch,
                                     const std::vector<int>& images) const {
  patch.m_images.clear();
  patch.m_grids.clear();
  vector<int>::const_iterator bimage = images.begin();
  vector<int>::const_iterator eimage = images.end();
  while (bimage != eimage) {
    const Vec3f icoord = m_fm.m_pss.project(*bimage, patch.m_coord, m_fm.m_level);
    const int ix = ((int)floor(icoord[0] + 0.5f)) / m_fm.m_csize;
    const int iy = ((int)floor(icoord[1] + 0.5f)) / m_fm.m_csize;
    if (0 <= ix && ix < m_gwidths[*bimage] &&
        0 <= iy && iy < m_gheights[*bimage]) {
      patch.m_images.push_back(*bimage);
      patch.m_grids.push_back(Vec2i(ix, iy));
    }
    ++bimage;
  }
}

void CpatchOrganizerS::setGrids(Ppatch& ppatch) const{
  setGrids(*ppatch);
}

void CpatchOrganizerS::setGrids(Cpatch& patch) const{
  patch.m_grids.clear();
  for (int i = 0; i < (int)patch.m_images.size(); ++i) {
    const int image = patch.m_images[i];
    Vec3f icoord = m_fm.m_pss.project(image, patch.m_coord, m_fm.m_level);
    const int ix = ((int)floor(icoord[0] + 0.5f)) / m_fm.m_csize;
    const int iy = ((int)floor(icoord[1] + 0.5f)) / m_fm.m_csize;
    patch.m_grids.push_back(TVec2<int>(ix, iy));
  }
}

void CpatchOrganizerS::setVImagesVGrids(Ppatch& ppatch) {
  setVImagesVGrids(*ppatch);
}

void CpatchOrganizerS::setVImagesVGrids(Cpatch& patch) {
  vector<int> used;
  used.resize(m_fm.m_tnum);
  fill(used.begin(), used.end(), 0);
  
  vector<int>::iterator bimage = patch.m_images.begin();
  vector<int>::iterator eimage = patch.m_images.end();
  while (bimage != eimage) {
    if ((*bimage) < m_fm.m_tnum)
      used[*(bimage)] = 1;
    ++bimage;
  }
  
  bimage = patch.m_vimages.begin();
  eimage = patch.m_vimages.end();
  while (bimage != eimage)
    used[*(bimage++)] = 1;

  for (int image = 0; image < m_fm.m_tnum; ++image) {
    if (used[image])
      continue;

    int ix, iy;
    if (isVisible0(patch, image, ix, iy,
                   m_fm.m_neighborThreshold, 1) == 0) {
      continue;
    }
    
    if (m_fm.m_pss.getEdge(patch.m_coord, image, m_fm.m_level) == 0)
      continue;
    
    patch.m_vimages.push_back(image);
    patch.m_vgrids.push_back(TVec2<int>(ix, iy));
  }
}

void CpatchOrganizerS::removePatch(const Ppatch& ppatch) {
  for (int i = 0; i < (int)ppatch->m_images.size(); ++i) {
    const int image = ppatch->m_images[i];
    if (m_fm.m_tnum <= image)
      continue;
    
    const int& ix = ppatch->m_grids[i][0];
    const int& iy = ppatch->m_grids[i][1];
    const int index = iy * m_gwidths[image] + ix;
    m_pgrids[image][index].erase(remove(m_pgrids[image][index].begin(),
                                        m_pgrids[image][index].end(),
                                        ppatch),
                                 m_pgrids[image][index].end());
  }

  for (int i = 0; i < (int)ppatch->m_vimages.size(); ++i) {
    const int image = ppatch->m_vimages[i];
#ifdef DEBUG
    if (m_fm.m_tnum <= image) {
      cerr << "Impossible in removePatch. m_vimages must be targetting images" << endl;
      exit (1);
    }
#endif

    const int& ix = ppatch->m_vgrids[i][0];
    const int& iy = ppatch->m_vgrids[i][1];
    const int index = iy * m_gwidths[image] + ix;
    m_vpgrids[image][index].erase(remove(m_vpgrids[image][index].begin(),
                                         m_vpgrids[image][index].end(),
                                         ppatch),
                                  m_vpgrids[image][index].end());
  }
}

int CpatchOrganizerS::isVisible0(const Cpatch& patch, const int image,
                                int& ix, int& iy,
                                const float strict, const int lock) {
  const Vec3f icoord =
    m_fm.m_pss.project(image, patch.m_coord, m_fm.m_level);
  ix = ((int)floor(icoord[0] + 0.5f)) / m_fm.m_csize;
  iy = ((int)floor(icoord[1] + 0.5f)) / m_fm.m_csize;

  return isVisible(patch, image, ix, iy, strict, lock);
}  

int CpatchOrganizerS::isVisible(const Cpatch& patch, const int image,
                               const int& ix, const int& iy,
                               const float strict, const int lock) {
  const int& gwidth = m_gwidths[image];
  const int& gheight = m_gheights[image];
  
  if (ix < 0 || gwidth <= ix || iy < 0 || gheight <= iy)
    return 0;

  if (m_fm.m_depth == 0)
    return 1;

  int ans = 0;
  Ppatch dppatch = m_MAXDEPTH;
  const int index = iy * gwidth + ix;
  
  if (lock)
    m_fm.m_imageLocks[image].rdlock();

  if (m_dpgrids[image][index] == m_MAXDEPTH)
    ans = 1;
  else
    dppatch = m_dpgrids[image][index];
  
  if (lock)
    m_fm.m_imageLocks[image].unlock();

  if (ans == 1)
    return 1;
  

  Vec4f ray = patch.m_coord - m_fm.m_pss.m_photos[image].m_center;
  unitize(ray);
  const float diff = ray * (patch.m_coord - dppatch->m_coord);
  const float factor = min(2.0, 2.0 + ray * patch.m_normal);
  
  if (diff < m_fm.m_optim.getUnit(image, patch.m_coord) * m_fm.m_csize * strict * factor)
    return 1;
  else
    return 0;
}  

void CpatchOrganizerS::findNeighbors(const Patch::Cpatch& patch,
                                     std::vector<Patch::Ppatch>& neighbors,
                                     const int lock, const float scale,
                                     const int margin,
                                     const int skipvis) {
  const float radius = 1.5 * margin * m_fm.m_expand.computeRadius(patch);
  
  vector<int>::const_iterator bimage = patch.m_images.begin();
  vector<int>::const_iterator eimage = patch.m_images.end();
  vector<Vec2i>::const_iterator bgrid = patch.m_grids.begin();

#ifdef DEBUG
  if (patch.m_images.empty()) {
    cerr << "Empty patches in findCloses" << endl;
    exit (1);
  }
#endif
  float unit = 0.0f;
  for (int i = 0; i < (int)patch.m_images.size(); ++i)
    unit += m_fm.m_optim.getUnit(patch.m_images[i], patch.m_coord);
  unit /= (int)patch.m_images.size();
  unit *= m_fm.m_csize;
  
  while (bimage != eimage) {
    if (m_fm.m_tnum <= *bimage) {
      ++bimage;
      ++bgrid;
      continue;
    }
    const int image = *bimage;
    const int& ix = (*bgrid)[0];
    const int& iy = (*bgrid)[1];
    if (lock)
      m_fm.m_imageLocks[image].rdlock();
    for (int j = -margin; j <= margin; ++j) {
      const int ytmp = iy + j;
      if (ytmp < 0 || m_fm.m_pos.m_gheights[image] <= ytmp)
        continue;
      for (int i = -margin; i <= margin; ++i) {
        const int xtmp = ix + i;
        if (xtmp < 0 || m_fm.m_pos.m_gwidths[image] <= xtmp)
          continue;
        const int index = ytmp * m_fm.m_pos.m_gwidths[image] + xtmp;
        vector<Ppatch>::const_iterator bpatch =
          m_fm.m_pos.m_pgrids[image][index].begin();
        vector<Ppatch>::const_iterator epatch =
          m_fm.m_pos.m_pgrids[image][index].end();
        while (bpatch != epatch) {
          if (m_fm.isNeighborRadius(patch, **bpatch, unit,
                                          m_fm.m_neighborThreshold * scale,
                                          radius))
            neighbors.push_back(*bpatch);
          ++bpatch;
        }
        bpatch = m_fm.m_pos.m_vpgrids[image][index].begin();
        epatch = m_fm.m_pos.m_vpgrids[image][index].end();
        while (bpatch != epatch) {
          if (m_fm.isNeighborRadius(patch, **bpatch, unit,
                                    m_fm.m_neighborThreshold * scale,
                                    radius))
            neighbors.push_back(*bpatch);
          ++bpatch;
        }
      }
    }
    if (lock)
      m_fm.m_imageLocks[image].unlock();

    ++bimage;
    ++bgrid;
  }

  if (skipvis == 0) {
    bimage = patch.m_vimages.begin();
    eimage = patch.m_vimages.end();
    bgrid = patch.m_vgrids.begin();
    
    while (bimage != eimage) {
      const int image = *bimage;
      const int& ix = (*bgrid)[0];
      const int& iy = (*bgrid)[1];
      if (lock)
        m_fm.m_imageLocks[image].rdlock();
      for (int j = -margin; j <= margin; ++j) {
        const int ytmp = iy + j;
        if (ytmp < 0 || m_fm.m_pos.m_gheights[image] <= ytmp)
          continue;
        for (int i = -margin; i <= margin; ++i) {
          const int xtmp = ix + i;
          if (xtmp < 0 || m_fm.m_pos.m_gwidths[image] <= xtmp)
            continue;
          const int index = ytmp * m_fm.m_pos.m_gwidths[image] + xtmp;
          vector<Ppatch>::const_iterator bpatch =
            m_fm.m_pos.m_pgrids[image][index].begin();
          vector<Ppatch>::const_iterator epatch =
            m_fm.m_pos.m_pgrids[image][index].end();
          while (bpatch != epatch) {
            if (m_fm.isNeighborRadius(patch, **bpatch, unit,
                                      m_fm.m_neighborThreshold * scale,
                                      radius))
              neighbors.push_back(*bpatch);
            ++bpatch;
          }
          bpatch = m_fm.m_pos.m_vpgrids[image][index].begin();
          epatch = m_fm.m_pos.m_vpgrids[image][index].end();
          while (bpatch != epatch) {
            if (m_fm.isNeighborRadius(patch, **bpatch, unit,
                                      m_fm.m_neighborThreshold * scale,
                                      radius))
              neighbors.push_back(*bpatch);
            ++bpatch;
          }
        }
      }
      if (lock)
        m_fm.m_imageLocks[image].unlock();
      
      ++bimage;
      ++bgrid;
    }
  }
  
  sort(neighbors.begin(), neighbors.end());
  neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
}

float CpatchOrganizerS::computeUnit(const Patch::Cpatch& patch) const{
  float unit = 0.0f;
  for (int i = 0; i < (int)patch.m_images.size(); ++i)
    unit += m_fm.m_optim.getUnit(patch.m_images[i], patch.m_coord);
  unit /= (int)patch.m_images.size();
  unit *= m_fm.m_csize;
  return unit;
}  

void CpatchOrganizerS::setScales(Patch::Cpatch& patch) const {
  const float unit = m_fm.m_optim.getUnit(patch.m_images[0], patch.m_coord);
  const float unit2 = 2.0f * unit;
  Vec4f ray = patch.m_coord - m_fm.m_pss.m_photos[patch.m_images[0]].m_center;
  unitize(ray);

  const int inum = min(m_fm.m_tau, (int)patch.m_images.size());
  
  // First compute, how many pixel difference per unit along vertical
  //for (int i = 1; i < (int)patch.m_images.size(); ++i) {
  for (int i = 1; i < inum; ++i) {
    Vec3f diff = m_fm.m_pss.project(patch.m_images[i], patch.m_coord, m_fm.m_level) -
      m_fm.m_pss.project(patch.m_images[i], patch.m_coord - ray * unit2, m_fm.m_level);
    patch.m_dscale += norm(diff);
  }

  // set m_dscale to the vertical distance where average pixel move is half pixel
  //patch.m_dscale /= (int)patch.m_images.size() - 1;
  patch.m_dscale /= inum - 1;
  patch.m_dscale = unit2 / patch.m_dscale;
  
  patch.m_ascale = atan(patch.m_dscale / (unit * m_fm.m_wsize / 2.0f));
}

void valid() {
  std::cout<<"valid"<<std::endl;
}

//Add Code
// algoritmo para transformar RGB em HSV
// parametro RGB na faixa de 0 a 255 e retorna HSV com H na faixa de 0 a 1.0 
/*vector<double>*/ Vec3f rgbToHsv(double r, double g, double b)
{
    r = r/255.0;
    g = g/255.0;
    b = b/255.0;

    double cmax = std::max(r, std::max(g, b));
    double cmin = std::min(r, std::min(g, b));
    double diff = cmax - cmin;
    double h = -1;
    double s = -1;

    if(cmax == cmin)
    {
        h = 0;
    }
    else if(cmax == r)
    {
        h = fmod(60.0 * ((g - b) / diff) + 360.0, 360.0);
    }
    else if(cmax == g)
    {
        h = fmod(60.0 * ((b - r) / diff) + 120.0, 360.0);
    }
    else if(cmax == b)
    {
        h = fmod(60.0 * ((r - g) / diff) + 240.0, 360.0);
    }

    if(cmax == 0)
    {
        s = 0;
    }
    else
    {
        s = (diff / cmax);
    }

    double v = cmax * 100;
     Vec3f colorf;
     colorf[0] = h /360;
     colorf[1] = s;
     colorf[2] = v;

    return colorf; //{h / 360.0, s, v};
}

// write out results
void CpatchOrganizerS::writePLY(const std::vector<Ppatch>& patches,
                                const std::string filename) {
  ofstream ofstr;
  ofstr.open(filename.c_str());
  ofstr << "ply" << '\n'
       << "format ascii 1.0" << '\n'
       << "element vertex " << (int)patches.size() << '\n'
       << "property float x" << '\n'
       << "property float y" << '\n'
       << "property float z" << '\n'
       << "property float nx" << '\n'
       << "property float ny" << '\n'
       << "property float nz" << '\n'
       << "property uchar diffuse_red" << '\n'
       << "property uchar diffuse_green" << '\n'
       << "property uchar diffuse_blue" << '\n'
       << "end_header" << '\n';

  vector<Ppatch>::const_iterator bpatch = patches.begin();
  vector<Ppatch>::const_iterator bend = patches.end();

  map<string,Vec3f> positionsColorsMap;

   int countImg = 0;
  while (bpatch != bend) {
    // Get color
    Vec3i color;
     countImg++;
    const int mode = 0;
    // 0: color from images
    // 1: fix
    // 2: angle
    if (mode == 0) {
      int denom = 0;
      int maxInterval = 4;
      int countBezier = 0;
      std::vector<Vec3f> colorBezier;

      Vec3f colorf;
      // if (countImg<300) {
      //   cout<<"######## before colors ##############"<<endl;
      // }

      int intervalMinMediaMovel = 4;
      int intervalMediaMovel = (int)(*bpatch)->m_images.size() < intervalMinMediaMovel? (int)(*bpatch)->m_images.size(): intervalMinMediaMovel;
      int countIntervalMediaMovel = 0;
      vector<Vec3f> valuesColorsMediaMovel;
      vector<Vec3f> valuesColorsResultMediaMovel;
      vector<Vec3f> valuesColorsList;
      vector<int> imageId;

      for (int i = 0; i < (int)(*bpatch)->m_images.size(); ++i) {
        const int image = (*bpatch)->m_images[i];
        imageId.push_back(image);//Ids das Imagens
        colorf += m_fm.m_Term_pss.getColor((*bpatch)->m_coord, image, m_fm.m_level);
        
        valuesColorsMediaMovel.push_back(m_fm.m_Term_pss.getColor((*bpatch)->m_coord, image, m_fm.m_level));

        valuesColorsList.push_back(m_fm.m_Term_pss.getColor((*bpatch)->m_coord, image, m_fm.m_level));

        if (countImg<300) {
           std::string vxx = std::to_string(12345);
           std::cout<< "init"<< vxx <<"########coordenada: "<< (*bpatch)->m_coord <<" image"<<image<<" level "<<m_fm.m_level<<std::endl;
        //   std::cout<<"######## color: "<<colorf<<" denom: "<<denom<<std::endl;
        }

        if (countIntervalMediaMovel < intervalMediaMovel) {
          countIntervalMediaMovel++;
        } else {
          int count = 0;
          Vec3f mediaMovel;
          while(count<intervalMediaMovel) {
            mediaMovel += valuesColorsMediaMovel[count];
            //std::cout<<"######## media movel: "<< mediaMovel[0] <<endl;
            count++;
          }
          //std::cout<<"######## media movel adicionada: "<< mediaMovel[0]<<endl;
          mediaMovel = mediaMovel/intervalMediaMovel;
          // std::cout<<"######## depois media movel adicionada: "<< mediaMovel[0]<<endl;
          valuesColorsResultMediaMovel.push_back(mediaMovel);
          
          valuesColorsMediaMovel.clear();          
          // std::cout<<"######## depois media movel adicionada: size "<< valuesColors.size() <<endl;
          countIntervalMediaMovel=0;
        }

        denom++;
      }

      //Experimento 1
      //colorf = valuesColorsResult[valuesColorsResult.size()-1];

      //Original
      // colorf /= denom;

      //Média movel
      Vec3f colorSumMediaMovel;
      for (int i=0; i< valuesColorsResultMediaMovel.size(); i++) {
        colorSumMediaMovel += valuesColorsResultMediaMovel[i];
      }
      colorf = colorSumMediaMovel/valuesColorsResultMediaMovel.size();

      //Bezier Início
      // Bezier::Bezier<4> poli({ {120, 160}, {35, 200}, {220, 260}, {220, 40}, { 210, 30} });
      
      //TODO: função para pegar tempo real e converter em minutos.
      int numImagens = (int)(*bpatch)->m_images.size();
      int timeSpendMinutes = numImagens * 5;
      float deltaT = timeSpendMinutes/timeSpendMinutes;

      vector<float> tPointList;
      for (int i=0; i< numImagens;i++) {
        tPointList.push_back(deltaT*imageId[i]);
      }

      vector<Vec3f> valuesColorsHsvList;

      for (int i=0; i< valuesColorsList.size(); i++) {
        Vec3f colorHsv = rgbToHsv(valuesColorsList[i][0], valuesColorsList[i][1], valuesColorsList[i][2]);
        valuesColorsHsvList.push_back(colorHsv);
      }      

      vector<Bezier::Point> hTvaluesBezierList;

      for (int i=0; i< numImagens; i++) {
        Vec2f valueHT;
        valueHT[0] = valuesColorsHsvList[i][0];
        valueHT[1] = tPointList[i];
        hTvaluesBezierList.push_back({ valuesColorsHsvList[i][0], tPointList[i] });
      }
      int qudList = hTvaluesBezierList.size();
      const int numPolinomios = 1>0? 2 : 3;  //(int)hTvaluesBezierList.size();
      Bezier::Bezier<qudList> poli({hTvaluesBezierList});

      for (float i = 0.0; i < 1.0f; i+=0.1) {
          Bezier::Point val = poli.valueAt(i);
          std::cout << " value P: " << i << " value 0  x " << val.x << " ..."<<std::endl;
          std::cout << " value 0 y " << val.y << " ..."<<std::endl;
      }

      float valueSelectedForPointT = 0.5; //TODO Selecionar via parametro
      Bezier::Point valuePointT = poli.valueAt(valueSelectedForPointT);

      // valuePointT.x

      //Bezier Fim
      
      color[0] = min(255,(int)floor(colorf[0] + 0.5f));
      color[1] = min(255,(int)floor(colorf[1] + 0.5f));
      color[2] = min(255,(int)floor(colorf[2] + 0.5f));

      // color[0] = 0.0f; //min(255,(int)floor(colorf[0] + 0.5f));
      // color[1] = 0.0f; //min(255,(int)floor(colorf[1] + 0.5f));
      // color[2] = 0.0f; //min(255,(int)floor(colorf[2] + 0.5f));

    }
    else if (mode == 1) {
      if ((*bpatch)->m_tmp == 1.0f) {
        color[0] = 255;
        color[1] = 0;
        color[2] = 0;
      }
      else {
        color[0] = 255;
        color[1] = 255;
        color[2] = 255;
      }
    }
    else if (mode == 2) {
      float angle = 0.0f;
      vector<int>::iterator bimage = (*bpatch)->m_images.begin();
      vector<int>::iterator eimage = (*bpatch)->m_images.end();

      while (bimage != eimage) {
        const int index = *bimage;
        Vec4f ray = m_fm.m_pss.m_photos[index].m_center - (*bpatch)->m_coord;
        ray[3] = 0.0f;
        unitize(ray);

        angle += acos(ray * (*bpatch)->m_normal);
        ++bimage;
      }
      
      angle = angle / (M_PI / 2.0f);
      float r, g, b;
      Image::Cimage::gray2rgb(angle, r, g, b);
      color[0] = (int)(r * 255.0f);
      color[1] = (int)(g * 255.0f);
      color[2] = (int)(b * 255.0f);
    }
    
    ofstr << (*bpatch)->m_coord[0] << ' '
          << (*bpatch)->m_coord[1] << ' '
          << (*bpatch)->m_coord[2] << ' '
          << (*bpatch)->m_normal[0] << ' '
          << (*bpatch)->m_normal[1] << ' '
          << (*bpatch)->m_normal[2] << ' '
          << color[0] << ' ' << color[1] << ' ' << color[2] << '\n';
      ++bpatch;
  }
  ofstr.close();  
}

void CpatchOrganizerS::writePLY(const std::vector<Ppatch>& patches,
                                const std::string filename,
                                const std::vector<Vec3i>& colors) {
  ofstream ofstr;
  ofstr.open(filename.c_str());
  ofstr << "ply" << '\n'
       << "format ascii 1.0" << '\n'
       << "element vertex " << (int)patches.size() << '\n'
       << "property float x" << '\n'
       << "property float y" << '\n'
       << "property float z" << '\n'
       << "property float nx" << '\n'
       << "property float ny" << '\n'
       << "property float nz" << '\n'
       << "property uchar diffuse_red" << '\n'
       << "property uchar diffuse_green" << '\n'
       << "property uchar diffuse_blue" << '\n'
       << "end_header" << '\n';

  vector<Ppatch>::const_iterator bpatch = patches.begin();
  vector<Ppatch>::const_iterator bend = patches.end();
  vector<Vec3i>::const_iterator colorb = colors.begin();
  
  while (bpatch != bend) {
    ofstr << (*bpatch)->m_coord[0] << ' '
          << (*bpatch)->m_coord[1] << ' '
          << (*bpatch)->m_coord[2] << ' '
          << (*bpatch)->m_normal[0] << ' '
          << (*bpatch)->m_normal[1] << ' '
          << (*bpatch)->m_normal[2] << ' '
          << *colorb << '\n';
    ++bpatch;
    ++colorb;
  }
  ofstr.close();  
}

//Add code
Vec3i HSVtoRGB(float H, float S,float V){
  Vec3i color;
  color[0] = 0;
  color[1] = 0;
  color[2] = 0;
    if(H>360 || H<0 || S>100 || S<0 || V>100 || V<0){
        cout<<"The givem HSV values are not in valid range"<<endl;
        return color;
    }
    float s = S/100;
    float v = V/100;
    float C = s*v;
    float X = C*(1-abs(fmod(H/60.0, 2)-1));
    float m = v-C;
    float r,g,b;
    if(H >= 0 && H < 60){
        r = C,g = X,b = 0;
    }
    else if(H >= 60 && H < 120){
        r = X,g = C,b = 0;
    }
    else if(H >= 120 && H < 180){
        r = 0,g = C,b = X;
    }
    else if(H >= 180 && H < 240){
        r = 0,g = X,b = C;
    }
    else if(H >= 240 && H < 300){
        r = X,g = 0,b = C;
    }
    else{
        r = C,g = 0,b = X;
    }
    int R = (r+m)*255;
    int G = (g+m)*255;
    int B = (b+m)*255;
    cout<<"R : "<<R<<endl;
    cout<<"G : "<<G<<endl;
    cout<<"B : "<<B<<endl;

    color[0] = R;
    color[0] = G;
    color[0] = B;

    return color; //{R, G, B};
}
