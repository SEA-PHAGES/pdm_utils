import os


def expand_path(input_path):
	"""
	Attempts to coerce input paths in any relative format into an
	absolute path.
	:param input_path: the path to be expanded
	:return expanded_path: the expanded path
	"""
	# "~" needs to be expanded first
	home_dir = os.path.expanduser("~")

	if input_path[0] == "~":
		expanded_path = home_dir + input_path[1:]
	else:
		expanded_path = input_path

	# os.path.abspath will resolve ./ or ../ present in the path
	expanded_path = os.path.abspath(expanded_path)
	return expanded_path


def verify_path(filepath, kind=None):
	"""
	Verifies that a given filepath is exists, and if a kind is given,
	it verifies that it exists as the indicated kind.
	:param filepath: full path to the desired file/directory.
	:param kind: ("file", "dir"), corresponding with paths to be
	checked as either files or directories.
	:return Boolean: True if filepath is verified, False otherwise.
	"""
	if kind == "file":
		if os.path.isfile(filepath) is True:
			return True
		else:
			return False
	elif kind == "dir":
		if os.path.isdir(filepath) is True:
			return True
		else:
			return False
	elif kind is None:
		if os.path.exists(filepath) is True:
			return True
		else:
			return False
	else:
		print("{} is not a valid kind for this function. Please try again "
			  "using one of (None, dir, file).")


def ask_yes_no(prompt):
	"""
	Function to get the user's yes/no response to a question.
	Accepts variations of yes/y, true/t, no/n, false/f.
	:param prompt: the question to ask the user.
	:return Boolean: default is False (e.g. user hits Enter w/o typing
	anything else), but variations of yes or true responses will return
	True instead.
	"""
	response = False
	response_valid = False
	while response_valid is False:
		response = input(prompt)
		if response.lower() in ["yes", "y", "t", "true"]:
			response = True
			response_valid = True
		elif response.lower() in ["no", "n", "f", "false", ""]:
			response = False
			response_valid = True
		else:
			print("Invalid response.")
	return response


def close_files(list_of_filehandles):
	"""
	Closes all the files in a list of open file handles.
	:param list_of_filehandles: A list of open file handles
	:return:
	"""
	for handle in list_of_filehandles:
		handle.close()
	return
