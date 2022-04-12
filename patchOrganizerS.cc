#define _USE_MATH_DEFINES
#include <cmath>

#include <string>
#include "patchOrganizerS.h"
#include "findMatch.h"

#include<bits/stdc++.h>
#include <set>
#include <vector>

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

Vec3f rgbToHsv(double r, double g, double b)
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

Vec3i hsvToRgb(float H, float S,float V) {
  Vec3i color;
  color[0] = 0;
  color[1] = 0;
  color[2] = 0;
    if(H>360 || H<0 || S>100 || S<0 || V>100 || V<0) {
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

double factorial(int num) {
    if (num <= 1) {
        return 1;
    } else {
        return num * factorial(num - 1);
    }
}

std::vector<double> bezier(double t, std::vector<std::vector<double> > bezierCtrlNodesArr) {

    double x = 0; 
    double y = 0;
    std::vector<std::vector<double> > bezierCtrlNodes = bezierCtrlNodesArr;

    int n = bezierCtrlNodesArr.size() - 1;
               
    for (int i = 0; i < bezierCtrlNodesArr.size(); i++) {
      try {
        std::cout << " x: " << bezierCtrlNodesArr[i][0] << " y: " << bezierCtrlNodesArr[i][1] << endl;
       if(!i) {
        x += bezierCtrlNodesArr[i][0] * pow((1 - t), n - i) * pow(t, i);
        y += bezierCtrlNodesArr[i][1] * pow((1 - t), n - i) * pow(t, i);
       } else {
        x += factorial(n) / factorial(i) / factorial(n - i) * bezierCtrlNodesArr[i][0] * pow(( 1 - t ), n - i) * pow(t, i);
        y += factorial(n) / factorial(i) / factorial(n - i) * bezierCtrlNodesArr[i][1] * pow(( 1 - t ), n - i) * pow(t, i);
       }
      } catch(...) {
        std::cout << "error i: " << i << " bezierCtrlNodesArr[i][0]: " << bezierCtrlNodesArr[i][0] << " bezierCtrlNodesArr[i][1]" << bezierCtrlNodesArr[i][1] << endl;
      }
    }

    std::vector<double> resultList;

    try{
    resultList.push_back(x);
    resultList.push_back(y);
    }catch(...) {
      std::cout << "error p2 "<< endl;
    }
    return resultList;
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

  //type get colors 0 -> media movel
  //type get colors 1 -> bezier
  //type get colors 2 -> media aritimética (padrão original)-

  //type get colors 3 -> primeira imagem do ponto (padrão)
  //type get colors 4 -> última imagem do ponto (padrão)
  //type get colors 5 -> imagem pelo tempo global (3 partes) 
  //type get colors 6 -> imagem pelo tempo global (2 partes)
  //type get colors 7 -> cada Imagem tem um tempo e cada ponto tem a média do tempo

  int typeGetColors = 7;
  vector <int> listImagemPreprcessing;
  vector <int> listImagem;
  int TR = 0; //Tempo de referência
  int timeImage = 1; //Tempo em minutos de cada imagem

  //mediaNeighborhood = 0 //Sem média na vizinhaça
  //mediaNeighborhood = 1 //Com média na vizinhaça
  int mediaNeighborhood = 0;
  
  vector<vector<double> > listValuesFile;

  vector<float> positionList;

  while (bpatch != bend) {
    positionList.push_back((*bpatch)->m_coord[2]);
    listImagemPreprcessing.push_back((*bpatch)->m_images[0]);
    ++bpatch;
  }

  //########################

  cout << "##################################listImagem:  " <<listImagem.size() << endl;

  std::set<int> list_set_preprocessing(listImagemPreprcessing.begin(), listImagemPreprcessing.end());
  listImagemPreprcessing.assign(list_set_preprocessing.begin(), list_set_preprocessing.end());
  
  cout << "list :  " <<listImagemPreprcessing.size() << endl;

  for (int i = 0; i < listImagemPreprcessing.size(); i++) {
    cout << "listImagemPreprcessing set list:  " <<listImagemPreprcessing[i] << endl;
  }
  cout << "listImagemPreprcessingset list size():  " <<listImagemPreprcessing.size() << endl;

  int totalImgPreprocessing = listImagemPreprcessing.size();
  
  TR = (totalImgPreprocessing * timeImage) / 2;

  //########################
  
  sort(positionList.begin(), positionList.end());

  int valorTotal = positionList.size();
  int tercaParte = valorTotal/3;

  int position = 0;
  int countList = 0;
  bpatch = patches.begin();

  int countImg = 0;
  while (bpatch != bend) {
    // Get color
    Vec3i color;
    countImg++;
    countList++;

    const int mode = 0;
    // 0: color from images
    // 1: fix
    // 2: angle
    if (mode == 0) {
      int denom = 0;
      Vec3f colorf;

      //int maxInterval = 4;
      int countBezier = 0;
      std::vector<Vec3f> colorBezier;

      int intervalMinMediaMovel = 4;
      int intervalMediaMovel = (int)(*bpatch)->m_images.size() < intervalMinMediaMovel? (int)(*bpatch)->m_images.size(): intervalMinMediaMovel;
      int countIntervalMediaMovel = 0;
      int numIntervalMediaMovel = 0;
      vector<Vec3f> valuesColorsResultMediaMovel;
      vector<Vec3f> valuesColorsList;
      vector<Vec3f> valuesMediaMovelColorsList;
      vector<int> imageId;

      vector<Vec3f> firstValue;
      vector<Vec3f> lastValue;

      vector<int> imageSortedList;

      for (int i = 0; i < (int)(*bpatch)->m_images.size(); ++i) {
        const int image = (*bpatch)->m_images[i];
        imageSortedList.push_back(image);        
      }
      sort(imageSortedList.begin(), imageSortedList.end());

      //cout << "lista ordenada Imagem size: " << imageSortedList.size() << endl;
      for (int i=0; i < imageSortedList.size(); i++) {

        if (i == 0) {
          //cout << " Primeiro valor: " << imageSortedList[i] << endl;
          firstValue.push_back(m_fm.m_Term_pss.getColor((*bpatch)->m_coord, imageSortedList[i], m_fm.m_level));
        }
        if (i+1 == imageSortedList.size()) {
          //cout << " Ultimo valor: " << imageSortedList[i] << endl;
          lastValue.push_back(m_fm.m_Term_pss.getColor((*bpatch)->m_coord, imageSortedList[i], m_fm.m_level));
        }
        Vec3f valuesColors = m_fm.m_Term_pss.getColor((*bpatch)->m_coord, imageSortedList[i], m_fm.m_level);
        valuesColorsList.push_back(valuesColors);        
      }

      for (int i = 0; i < (int)(*bpatch)->m_images.size(); ++i) {
        const int image = (*bpatch)->m_images[i];
        imageId.push_back(image);
        colorf += m_fm.m_Term_pss.getColor((*bpatch)->m_coord, image, m_fm.m_level);

        Vec3f valuesColors = m_fm.m_Term_pss.getColor((*bpatch)->m_coord, image, m_fm.m_level);

        valuesMediaMovelColorsList.push_back(valuesColors);

        if (typeGetColors == 0) {
          if (countIntervalMediaMovel < (intervalMediaMovel -1)) {
            std::cout << " countIntervalMediaMovel: "<< countIntervalMediaMovel << "intervalMediaMovel: " << intervalMediaMovel << endl;
            countIntervalMediaMovel++;
          } else {

            std::cout << "CountInterval MediaMovel: "<< countIntervalMediaMovel << " IMAGEM SIZE: " << (int)(*bpatch)->m_images.size() << endl;
            Vec3f mediaMovel;
            mediaMovel = 0;
            int intervalInit = countIntervalMediaMovel - (intervalMediaMovel-1);
            int intervalFinish = intervalInit + intervalMediaMovel;

            std::cout << " intervalInit: " << intervalInit << " finish: " << intervalFinish << "intervalMediaMovel: " << intervalMediaMovel << endl;

            for(int index = intervalInit; index < intervalFinish; index++) {
              mediaMovel += valuesMediaMovelColorsList[index];
            }
            mediaMovel = mediaMovel/intervalMediaMovel;
            std::cout << " mediaMovel: " << mediaMovel << endl;
            valuesColorsResultMediaMovel.push_back(mediaMovel);
            
            countIntervalMediaMovel++;
            
          }
        }

        denom++;
      }

      if (typeGetColors == 0) {
        Vec3f colorSumMediaMovel;
        for (int i=0; i< valuesColorsResultMediaMovel.size(); i++) {
          colorSumMediaMovel += valuesColorsResultMediaMovel[i];
        }
        //Passa valor da média móvel para o final
        colorf = colorSumMediaMovel/valuesColorsResultMediaMovel.size();
      }

      if (typeGetColors == 1) {
        cout << " bezier" << endl;
        int numImagens = (int)(*bpatch)->m_images.size();
        //int timeSpendMinutes = numImagens * 5 * 60;
        int timeSpendMinutes = 60;
        float deltaT = timeSpendMinutes/numImagens;

        vector<float> tPointList;
        for (int i=0; i< numImagens;i++) {
          tPointList.push_back(deltaT*imageId[i]);
        }
        sort(tPointList.begin(), tPointList.end());

        vector<Vec3f> valuesColorsHsvList;
        for (int i=0; i< valuesColorsList.size(); i++) {
          valuesColorsHsvList.push_back(rgbToHsv(valuesColorsList[i][0], valuesColorsList[i][1], valuesColorsList[i][2]));
        }

        std::vector<std::vector<double> > valuesR;
        for (int i=0; i<numImagens; i++) {
          std::vector<double> valueHT;
          std::vector<double> valueHTHsv;
          setprecision(2);

          valueHT.push_back(round(valuesColorsList[i][0]));
          valueHTHsv.push_back(round(valuesColorsHsvList[i][0]));
          valueHT.push_back(tPointList[i]);
          valuesR.push_back(valueHT);
          // valuesR.push_back(valueHTHsv);
        }

        std::vector<std::vector<double> > valuesG;
        for (int i=0; i<numImagens; i++) {
          std::vector<double> valueHT;
          std::vector<double> valueHTHsv;
          setprecision(2);

          valueHT.push_back(round(valuesColorsList[i][1]));
          valueHTHsv.push_back(round(valuesColorsHsvList[i][0]));
          valueHT.push_back(tPointList[i]);
          valuesG.push_back(valueHT);
        }

        std::vector<std::vector<double> > valuesB;
        for (int i=0; i<numImagens; i++) {
          std::vector<double> valueHT;
          setprecision(2);

          valueHT.push_back(round(valuesColorsList[i][2]));
          valueHT.push_back(tPointList[i]);
          valuesB.push_back(valueHT);
        }            

        double calc = 0.0;
        for(int i = 0; i < tPointList.size(); i++) {
          cout << "tPointList[i]: "<< tPointList[i] << " R: " << valuesColorsList[i][0] << " G: " << valuesColorsList[i][1] << " B: " << valuesColorsList[i][2] << endl;
          double ai = (tPointList[i] - tPointList[0]) / (tPointList[tPointList.size()-1] - tPointList[0]);
          cout << " ai " << ai << endl;
          calc += ai;
        }

        double paramT = calc/tPointList.size();//0.4f;
        setprecision(2);

        cout << " calc " << calc << " paramT " << paramT << " tPointList.size() " << tPointList.size() << endl;
        
        std::vector<double> resultBezierR;
        std::vector<double> resultBezierG;
        std::vector<double> resultBezierB;

        //if (valuesR.size() <31) {
          resultBezierR = bezier(paramT, valuesR);
          cout << " G: " << endl;
          resultBezierG = bezier(paramT, valuesG);
          cout << " B: " << endl;
          resultBezierB = bezier(paramT, valuesB);

        cout<< "paramT:"<< paramT << endl;

        colorf[0] = resultBezierR[0];
        colorf[1] = resultBezierG[0];
        colorf[2] = resultBezierB[0];              
      }

      if (typeGetColors == 2) {
        cout << " média aritimética" << endl;
        colorf /= denom;
      }

      if (typeGetColors == 3) {
        cout << "Primeiro valor, size: " << firstValue.size() << endl;
        colorf = firstValue[0];

      }

      if (typeGetColors == 4) {
        
        cout << "C++:" << __cplusplus << " Último valor, size: "<< lastValue.size() << endl;
        colorf = lastValue[0];
      }      

      if (typeGetColors == 5) {
        
        int imageProcessing = (*bpatch)->m_images[0];
        int mediaImagem = 0;

        for (int i = 0; i < (*bpatch)->m_images.size(); i++) {
          cout << "imagen["<< i <<"]: "<< (*bpatch)->m_images[i] << endl;
          mediaImagem += (*bpatch)->m_images[i];
        }

        mediaImagem = mediaImagem/(*bpatch)->m_images.size();

        cout << "Num Imagem: "<< (*bpatch)->m_images.size() << endl;
        cout << "Tempo global: "<< valuesColorsList.size() << endl;
        cout << "Total imagem: " << totalImgPreprocessing << endl;
        listImagem.push_back((*bpatch)->m_images[0]);

        //TODO: 
        // (1) Dividir em dois pedaços
        // (2) Dividir em 3 (inverter o que tem hoje, primeiros pegar os últimas imagens)
        // (3) Pegar a imagem mais próxima ao tempo de referência, por parâmetro(pode ser uma porcentagem, podendo ser 0 inicial e 1 tempo final)
        // (4) Histograma

        if (mediaImagem <= (0.4 * totalImgPreprocessing)) {
          colorf = lastValue[0];
          //colorf = firstValue[0]; //valuesColorsList[0];
          cout << "0.4 mediaImagem:"<< mediaImagem << " Imagem:" << (*bpatch)->m_images[0] << endl;
        } else if (mediaImagem > (0.4 * totalImgPreprocessing) 
        && mediaImagem < (0.6 * totalImgPreprocessing)) {
          int tamanhoListaCores = valuesColorsList.size();
          int posicao = tamanhoListaCores / 2;
          colorf = valuesColorsList[posicao];
          cout << "0.4> x <0.6 mediaImagem:"<< mediaImagem << " Imagem:" << (*bpatch)->m_images[0] << endl;
        } else if (mediaImagem >= (0.6 * totalImgPreprocessing)){
          colorf = firstValue[0];
          //colorf = lastValue[0]; //valuesColorsList[valuesColorsList.size() -1];
          cout << "0.6 mediaImagem:"<< mediaImagem << " Imagem:" << (*bpatch)->m_images[0] << endl;
        }
      }

      if (typeGetColors == 6) {

        int mediaImagem = 0;
        // for (int i = 0; i < (*bpatch)->m_images.size(); i++) {
        //   cout << "imagen["<< i <<"]: "<< (*bpatch)->m_images[i] << endl;
        //   mediaImagem += (*bpatch)->m_images[i];
        // }

        //mediaImagem = mediaImagem/(*bpatch)->m_images.size();

        cout << "Num Imagem: "<< (*bpatch)->m_images.size() << endl;
        cout << "Tempo global: "<< valuesColorsList.size() << endl;
        cout << "Total imagem: " << totalImgPreprocessing << endl;
        listImagem.push_back((*bpatch)->m_images[0]);

        //TODO: 
        // (1) Dividir em dois pedaços
        // (2) Dividir em 3 (inverter o que tem hoje, primeiros pegar os últimas imagens)
        // (3) Pegar a imagem mais próxima ao tempo de referência, por parâmetro(pode ser uma porcentagem, podendo ser 0 inicial e 1 tempo final)
        // (4) Histograma

        if ((*bpatch)->m_images[0] <= (0.5 * totalImgPreprocessing)) {
          //colorf = firstValue[0]; //valuesColorsList[0];
          colorf = lastValue[0];
          cout << "menor 0.5: "<< (*bpatch)->m_images[0] << endl;
        } else { //if ((totalImgPreprocessing - mediaImagem) > (0.5 * totalImgPreprocessing)){
          //colorf = lastValue[0]; //(*bpatch)->m_images[(*bpatch)->m_images.size()-1]; //valuesColorsList[valuesColorsList.size() -1];
          colorf = firstValue[0];
          cout << "maior 0.5: "<< (*bpatch)->m_images[0] << endl;
        }
      } else if (typeGetColors == 7) {

        int timeAvaregeImage = 0;
        int timeImageValue = 0;
        Vec3f imageNear;
        int nearTR = -1;
        int diff = 0;
        timeImage = 1;

        for (int i = 0; i < (*bpatch)->m_images.size(); i++) {
          
          timeImageValue = (*bpatch)->m_images[i] * timeImage;
          timeAvaregeImage += timeImageValue;
          cout << "imagen["<< i <<"]: "<< (*bpatch)->m_images[i] << " timeImageValue "<< timeImageValue <<" timeImage: " << timeImage << "timeAvaregeImage:" << timeAvaregeImage << endl;
          // diff = (timeImageValue - TR);
          // diff = diff < 0? diff * -1: diff;
          // if (nearTR == -1) {
          //   nearTR = diff;
          // }
          // if (nearTR > diff) {
          //   nearTR = diff;
          //   imageNear = m_fm.m_Term_pss.getColor((*bpatch)->m_coord, (*bpatch)->m_images[i], m_fm.m_level);
          // }
        }

        //colorf = imageNear;
        
        //Média da imagem
        timeAvaregeImage = timeAvaregeImage/(*bpatch)->m_images.size();

        cout << "TR " << TR << "timeAvaregeImage" << timeAvaregeImage << endl;

        if (timeAvaregeImage < TR) {
          int imagem = (int)((*bpatch)->m_images[(*bpatch)->m_images.size()-1] * timeImage);
          if (imagem < TR) {
            colorf = lastValue[0];
          } else {
            timeImageValue = 0;
            diff = 0;
            nearTR = -1;
            for (int i = 0; i < (*bpatch)->m_images.size(); i++) {
              timeImageValue = (*bpatch)->m_images[i] * timeImage;
              diff = (timeImageValue - TR);
              diff = diff < 0? diff * -1: diff;
              if (nearTR == -1) {
                nearTR = diff;
                imageNear = lastValue[0];
              } 
              if (nearTR > diff) {
                nearTR = diff;
                imageNear = m_fm.m_Term_pss.getColor((*bpatch)->m_coord, (*bpatch)->m_images[i], m_fm.m_level);
              }
            }
            colorf = imageNear;
          }
            
        } else {
          int imagem = (int)((*bpatch)->m_images[0] * timeImage);

          if (imagem > TR) {
            colorf = firstValue[0];
          } else {
            timeImageValue = 0;
            timeImage = 0;
            diff = 0;
            nearTR = -1;
            for (int i = 0; i < (*bpatch)->m_images.size(); i++) {
              timeImageValue = (*bpatch)->m_images[i] * timeImage;
              diff = (timeImageValue - TR);
              diff = diff < 0? diff * -1: diff;
              //cout << "diff : "<< diff << " imagem: " << (*bpatch)->m_images[i] << endl;
              if (nearTR == -1) {
                nearTR = diff;
                imageNear = firstValue[0];
              } 
              if (nearTR > diff) {
                //cout << " nearTR: " << nearTR << "diff : "<< diff << endl;
                nearTR = diff;
                imageNear = m_fm.m_Term_pss.getColor((*bpatch)->m_coord, (*bpatch)->m_images[i], m_fm.m_level);
              }
            }
            colorf = imageNear;
          }
        }
      }

      color[0] = min(255,(int)floor(colorf[0] + 0.5f));
      color[1] = min(255,(int)floor(colorf[1] + 0.5f));
      color[2] = min(255,(int)floor(colorf[2] + 0.5f));
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
      float angle = 0.5f;
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

    vector<double> listValuesRow;

    listValuesRow.push_back((*bpatch)->m_coord[0]);
    listValuesRow.push_back((*bpatch)->m_coord[1]);
    listValuesRow.push_back((*bpatch)->m_coord[2]);

    listValuesRow.push_back((*bpatch)->m_normal[0]);
    listValuesRow.push_back((*bpatch)->m_normal[1]);
    listValuesRow.push_back((*bpatch)->m_normal[2]);

    listValuesRow.push_back(color[0]);
    listValuesRow.push_back(color[1]);
    listValuesRow.push_back(color[2]);

    listValuesFile.push_back(listValuesRow);
    
    if (mediaNeighborhood == 0) {
      ofstr << (*bpatch)->m_coord[0] << ' '
            << (*bpatch)->m_coord[1] << ' '
            << (*bpatch)->m_coord[2] << ' '
            << (*bpatch)->m_normal[0] << ' '
            << (*bpatch)->m_normal[1] << ' '
            << (*bpatch)->m_normal[2] << ' '
            << color[0] << ' ' << color[1] << ' ' << color[2] << '\n';
    }

      ++bpatch;
  }

  cout << "##################################listImagem:  " <<listImagem.size() << endl;

  std::set<int> list_set(listImagem.begin(), listImagem.end());
  listImagem.assign(list_set.begin(), list_set.end());
  
  cout << "list :  " <<listImagem.size() << endl;

  for (int i = 0; i< listImagem.size(); i++) {
  cout << "set list:  " <<listImagem[i] << endl;
  }
  cout << "set list size():  " <<listImagem.size() << endl;


  if (mediaNeighborhood == 1) {
    // int countGroup = 0;
    // int countGroupNF = 0;
    //Lista //Grupo //Linha
    vector<vector<vector<double> > > groupList;
    for(int i = 0; i < listValuesFile.size(); i++) {
      if (groupList.size()==0) {
        //Grupo//Linha
        vector<vector<double> > group;
        group.push_back(listValuesFile[i]);
        groupList.push_back(group);

      } else {
        
        bool addGroup = false;
        for (int j=0; j < groupList.size(); j++) {
          
          double space = groupList[j][0][0] - listValuesFile[i][0];
          space = space < 0? space * -1: space;

          if (space < 0.0003) {
           groupList[j].push_back(listValuesFile[i]);
            addGroup= true;
          }

        }

        if (!addGroup) {
          // countGroupNF++;
          // if (countGroupNF< 300) {
          //   cout << "#############Grupo não encontrado cria: " << endl;
          // }
          
          vector<vector<double> > groupedValues;
          groupedValues.push_back(listValuesFile[i]);
          groupList.push_back(groupedValues);        
        }
      }
    }

    cout << "1 groupList: " << groupList.size() << " first group rows: " << groupList[0].size() << " col group: " << groupList[0][0].size() << endl;

    int intervalMediaMovel = 0;
    int maxIntervalMediaMovel = 4;
    vector <double> mediaMovelList;
    vector <double> mediaMovelGList;
    vector <double> mediaMovelBList;
    vector <double> mediaMoveResultlList;
    vector <double> mediaMoveResultlGList;
    vector <double> mediaMoveResultlBList;
    
    int countValues = 0;

    for(int i = 0; i < groupList.size(); i++) {
      for ( int j = 0; j < groupList[i].size(); j++) {

        mediaMovelList.push_back(groupList[i][j][6]);
        mediaMovelGList.push_back(groupList[i][j][7]);
        mediaMovelBList.push_back(groupList[i][j][8]);                  

        if (intervalMediaMovel < (maxIntervalMediaMovel-1)) {
          countValues++;

          if (countValues< 30) {
            cout << " Adiciona: " << intervalMediaMovel << " maxIntervalMediaMovel+intervalMediaMovel: " << maxIntervalMediaMovel+intervalMediaMovel << " groupList[i][j][6] "<< groupList[i][j][6]<< " groupList[i][j][7] " << groupList[i][j][7] << endl;
          }
          
          //Get R

          intervalMediaMovel +=1;
        } else {
          if (countValues< 30) {
            cout << "..else media" << endl;
          }
          
          int sumValuesR = 0; 
          int sumValuesG = 0; 
          int sumValuesB = 0; 

          int initInterval = intervalMediaMovel - (maxIntervalMediaMovel-1);
          int finishInterval = initInterval + maxIntervalMediaMovel;

          if (countValues< 30) {
            cout << " initInterval: " << initInterval << " finishInterval: " << finishInterval << endl;
          }

          for (int k = initInterval; k < finishInterval; k++) {

            if (countValues< 30) {
              cout << " mediaMovelListR: " << mediaMovelList[k] << " G: "<< mediaMovelGList[k]<< "B: "<< mediaMovelBList[k] << endl;
            }            

            sumValuesR += mediaMovelList[k];
            sumValuesG += mediaMovelGList[k];
            sumValuesB += mediaMovelBList[k];
          }

          int mediaR = sumValuesR/maxIntervalMediaMovel;
          int mediaG = sumValuesG/maxIntervalMediaMovel;
          int mediaB = sumValuesB/maxIntervalMediaMovel;

          if (countValues< 30) {
            cout << " mediaR: " << mediaR << " mediaG: " << mediaG << " mediaB " << mediaB << endl;
          }
          
          
          mediaMoveResultlList.push_back(mediaR);
          mediaMoveResultlGList.push_back(mediaG);
          mediaMoveResultlBList.push_back(mediaB);
          intervalMediaMovel +=1;
        }
      }

      int mediaR = 0;
      int mediaG = 0;
      int mediaB = 0;
      for(int l=0; l< mediaMoveResultlList.size();l++) {
        mediaR += mediaMoveResultlList[i];
        mediaG += mediaMoveResultlGList[i];
        mediaB += mediaMoveResultlBList[i];
      }
      mediaR = mediaR/mediaMoveResultlList.size();
      mediaG = mediaG/mediaMoveResultlList.size();
      mediaB = mediaB/mediaMoveResultlList.size();

      if (countValues< 30) {
        cout << " mediaR: " << mediaR << " mediaG: " << mediaG << " mediaB " << mediaB << endl;
      }

      for ( int j = 0; j < groupList[i].size(); j++) {
        groupList[i][j][6] = mediaR;
        groupList[i][j][7] = mediaG;
        groupList[i][j][8] = mediaB;
      }
    }

    //Fim Média
    cout << "1##################################listImagem:  " <<listImagem.size() << endl;
    //Carrega   valores
    int countVl = 0;
    for(int i = 0; i < listValuesFile.size(); i++) {
      for(int j = 0; j < groupList.size(); j++) {
        bool finded = false;
        for ( int k = 0; k < groupList[j].size(); k++) {
          
          // if (countVl < 30) {
          //     cout << "listValuesFile[i][0] " << listValuesFile[i][0] << " groupList[j][k][0]: " << groupList[j][k][0] << endl;
          // }

          if (listValuesFile[i][0] == groupList[j][k][0]) {
            countVl ++;
            if (countVl < 30) {
              cout << "finded: listValuesFile[i][0] " << listValuesFile[i][0] << " groupList[j][k][0] " << groupList[j][k][0] << (listValuesFile[i][0] == groupList[j][k][0]) << endl;
              cout << "groupList[j][k][6]:" << groupList[j][k][6] << " groupList[j][k][7] " << groupList[j][k][7] << " groupList[j][k][8]" << groupList[j][k][8] << endl;
            }
            listValuesFile[i][6] = groupList[j][k][6];
            listValuesFile[i][7] = groupList[j][k][7];
            listValuesFile[i][8] = groupList[j][k][8];
            finded = true;
            break;
          }
        }
        if (finded){
          break;
        }
      }
    }

  cout << "##################################listImagem:  " <<listImagem.size() << endl;

  std::set<int> list_set(listImagem.begin(), listImagem.end());
  listImagem.assign(list_set.begin(), list_set.end());
  
  cout << "list :  " <<listImagem.size() << endl;

  for (int i = 0; i< listImagem.size(); i++) {
    cout << "set list:  " <<listImagem[i] << endl;
  }
  cout << "...set list size():  " <<listImagem.size() << endl;
   

    cout << "write file" << endl;
    //Escreve no arquivo
    for(int i = 0; i < listValuesFile.size(); i++) {
      ofstr << listValuesFile[i][0] << ' '
            << listValuesFile[i][1] << ' '
            << listValuesFile[i][2] << ' '
            << listValuesFile[i][3] << ' '
            << listValuesFile[i][4] << ' '
            << listValuesFile[i][5] << ' '
            << listValuesFile[i][6] << ' ' << listValuesFile[i][7] << ' ' << listValuesFile[i][8] << '\n';

      // ofstr << (*bpatch)->m_coord[0] << ' '
      //       << (*bpatch)->m_coord[1] << ' '
      //       << (*bpatch)->m_coord[2] << ' '
      //       << (*bpatch)->m_normal[0] << ' '
      //       << (*bpatch)->m_normal[1] << ' '
      //       << (*bpatch)->m_normal[2] << ' '
      //       << color[0] << ' ' << color[1] << ' ' << color[2] << '\n';

    }
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
