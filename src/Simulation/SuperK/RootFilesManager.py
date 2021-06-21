import ROOT as root

def read_tree(filename, treename):
	file = root.TFile.Open(filename,'read')
	tree = file.Get(treename)
	#tree.Print()
	return file, tree

def out_tree(filename, treename):
	file = root.TFile.Open(filename,'recreate')
	tree = root.TTree(treename,treename)
	#tree.Print()
	return file, tree

def read_histo(filename):
	file = root.TFile.Open(filename,'read')
	return file

# def out_tree_variables(tree, variable_name, dimension):
# 	var  = np.zeros(dimension,'double')

# 	tree.Branch(variable_name,var,variable_name+'/D')

# 	return tree



# def root_file_histo:
# 	def __init__(self, filename):
# 		self.filename = filename	