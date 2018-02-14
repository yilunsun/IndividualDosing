#ifndef DELTA_H
#define DELTA_H

#ifndef _DIAGRAMS_DELTA_H
#define _DIAGRAMS_DELTA_H

#define DELTA_BASIC  1
#define DELTA_CHANGE 2
#define DELTA_PRUNE  3
#define DELTA_GROW   4

class Delta 
{
 public:
  Delta(int type);
  virtual ~Delta(void);

  virtual NodeTree *PerformChange(NodeTree *tree) const = 0;
  int GetType(void) const {return type;};
  
 private:
  int type;
};


inline Delta::Delta(int type)
{
  this->type = type;

  return;
};


inline Delta::~Delta(void)
{
  return;
};


#endif
#endif