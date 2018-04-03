def GetNameWOChain(whole_name):#canviar nom amb SUbl
	""" 
	Gets the name without the chain of a protein ID (name will be of the format "abc").
	
	Arguments:

	whole_name: name of the protein in format "abc_A"
	"""

	return whole_name[:-2]

def GetTargInteractionKeys(Final_interactions):
	"""
	Gets the keys of target interactions (from Final_interactions dictionary).
	"""

	return Final_interactions["target_interacts"].keys()

def GetTempInteractionKeys(Final_interactions, temp):
	"""
	Gets the keys of template interactions (from Final_interactions dictionary).
	"""
	return Final_interactions["temps"][temp]["temp_interact"].keys()

def GetChain(whole_name):
	""" 
	Gets the chain name of a protein ID (chain name will be of the format "A").
	
	Arguments:

	whole_name: name of the protein of format "abc_A"
	"""

	return whole_name[-1:]

def GetTempInteractions(Final_interactions, temp_name, temp):
	"""
	Gets the values of the template interactions (from Final_interactions dictionary).
	"""

	return Final_interactions["temps"][temp]["temp_interact"][GetChain(temp_name)]

def GetTargetInteractions(Final_interactions, targ_name):
	"""
	Gets the values of the target interactions (from Final_interactions dictionary).
	"""

	return Final_interactions["target_interacts"][GetChain(targ_name)]
