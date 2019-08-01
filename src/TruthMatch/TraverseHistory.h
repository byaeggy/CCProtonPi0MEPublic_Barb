#ifndef TraverseHistory_h_seen
#define TraverseHistory_h_seen

#include <map>
#include <vector>

namespace Minerva {
    class TG4Trajectory;
}

using Minerva::TG4Trajectory;

class TraverseHistory {
  public:
    explicit TraverseHistory(int trackId);

    void DoTraverse(const std::map<int, TG4Trajectory*>& trajectoryMap);
    void DoTraverse(const std::map<int, TG4Trajectory*>& trajectoryMap, int goalcode);

    int GetPDGCode() const;
    int GetPrimaryId() const;
    int GetPrimaryPdg() const;
    int GetParentId() const;
    int GetParentPdg() const;
    int GetNextToPrimaryId() const;
    int GetNextToPrimaryPdg() const;
    bool IsPresent(int pdg) const;
    
  private:
    int fTrackId;
    int fPdg;
    std::vector<int> fHistoryTrackId;
    std::vector<int> fHistoryPdg;
};


#endif
