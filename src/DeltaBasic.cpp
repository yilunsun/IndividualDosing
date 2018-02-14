#include "Random.h"
#include "Observation.h"
#include "Density.h"
#include "Model.h"
#include "NodeTree.h"
#include "DeltaBasic.h"


DeltaBasic::DeltaBasic(NodeTree *tree) : Delta(DELTA_BASIC)
{
	newTree = tree;

	return;
}



DeltaBasic::~DeltaBasic(void)
{
	newTree->RemoveSubTree();

	delete newTree;

	return;
};



NodeTree *DeltaBasic::PerformChange(NodeTree *tree) const
{
	tree->RemoveSubTree();
	delete tree;
	tree = NULL;
	return newTree->CopyTree();
};
