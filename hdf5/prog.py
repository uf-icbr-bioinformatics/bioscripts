#!/usr/bin/env python3

import sys
import os
import regex
import h5py
import math

class SeqBuilder:
	def __init__(self, seq_file_name, token_to_search, hdf5_file_name):
		self.seq_file_name = seq_file_name
		self.token_to_search = token_to_search
		self.hdf5_file_name = hdf5_file_name
		self.token_expr = regex.compile(token_to_search)
		self.max_token_length = 3

	#*******************************************************************************
	#*******************************************************************************

	def get_header_line (self, seq_file, hdf5_file):
		chrom_name = ""
		header_line_char = None
		while ((header_line_char := seq_file.read (1)) != ''):

			#compose header line
			chrom_name = chrom_name + header_line_char

			#when end of header line is reached,
			#note down the next character position as sequence starting position
			if (header_line_char == "\n"):
				seq_begin = seq_file.tell () + 1
				chrom_name = chrom_name.strip ()
				hdf5_file.create_group(chrom_name)
				return (chrom_name, seq_begin)

	#*******************************************************************************
	#*******************************************************************************

	def get_bucket_name (self, curr_position):
		return str (math.ceil (curr_position / 5000) * 5000)

	#*******************************************************************************
	#*******************************************************************************

	def add_token_to_hdf5_file (self, hdf5_file, chrom_name, curr_token, curr_position):
		curr_keys = hdf5_file[chrom_name].keys()

		#If the current token (say "CG") group is not yet in created Hd5File , create it
		if (curr_token not in list(curr_keys)):
			hdf5_file[chrom_name].create_group (curr_token)

		#for performance reasons - divide every 5000 indexes into a separate bucket
		#0-5000: will be stored in CHROMOSOME_NAME/SEQ_NAME/5000
		#5001-10000: will be stored in CHROMOSOME_NAME/SEQ_NAME/10000   etc.
		bucket_num_str = self.get_bucket_name (curr_position)
		if (bucket_num_str not in list (hdf5_file[chrom_name+"/"+curr_token].keys ())):
			hdf5_file.create_dataset (chrom_name+"/"+curr_token+"/"+bucket_num_str, (5000,),
									dtype='i4', compression="lzf", chunks=1000, fillvalue=-1)

		curr_dataset = hdf5_file[chrom_name+"/"+curr_token+"/"+bucket_num_str]

		if(not curr_dataset.attrs.__contains__('max_index')):
			curr_dataset.attrs['max_index'] = 0

		#if (curr_dataset.attrs['max_index'] % 20 == 0):
		print (chrom_name, curr_token, curr_position)

		curr_dataset[curr_dataset.attrs['max_index']] = curr_position

		curr_dataset.attrs['max_index'] += 1

	#*******************************************************************************
	#*******************************************************************************

	def process_tokens (self):
		new_line = True
		curr_token = ""
		curr_char = None
		chrom_name = None
		seq_begin = -1
		curr_position = 0
		match_obj = None
		
		#open the chromosomes/sequence file
		curr_file = open(self.seq_file_name, 'r')
		hdf5_file = h5py.File (self.hdf5_file_name, "w")
		
		while ((curr_char := curr_file.read (1)) != ''):

			#ignore new line characters
			if (curr_char == "\n"):
				new_line = True
				continue

			#If you got to a new Chromosome - return to main function
			if (new_line and curr_char == ">"):
				(chrom_name, seq_begin) = self.get_header_line (curr_file, hdf5_file)
				continue

			curr_token += curr_char.upper ()

			while (len(curr_token) > 0):
				#If we found a match for a full sequence - add to file
				if((match_obj := self.token_expr.fullmatch (curr_token)) != None):
					self.add_token_to_hdf5_file (hdf5_file, chrom_name, curr_token, curr_position-seq_begin)

					#If this sequence is part of a GROUP (like CHG for CAG) - add group to file
					grp_dict = match_obj.groupdict () 
					for gg in grp_dict.keys ():
						if (grp_dict[gg] == curr_token):
							self.add_token_to_hdf5_file (hdf5_file, chrom_name, gg, curr_position-seq_begin)

				#if there is no partial match for current token - try right part of token
				#CTG also has TG as part of that token
				if (not self.token_expr.fullmatch (curr_token, partial=True)):
					curr_token = curr_token[1:]
					curr_position += 1
				else:
					#if the first letter of a sequence is found - note down the position
					if (len (curr_token) == 1):
						curr_position = curr_file.tell ()
					break

		curr_file.close ()
		hdf5_file.close ()

	#*******************************************************************************
	#*******************************************************************************

	def check_token (self, chromosome_name, token_to_search, index_num):

		hdf5_file = h5py.File (self.hdf5_file_name, "r")
		if (chromosome_name not in list(hdf5_file.keys())):
			print ("Did not find Chromosome: " + chromosome_name)
			return
		else:
			print ("Found Chromosome: " + chromosome_name)
		
		if (token_to_search not in list(hdf5_file[chromosome_name].keys())):
			print ("Did not find Sequence: " + token_to_search)
			return
		else:
			print ("Found Sequence: " + token_to_search)

		bucket_name = self.get_bucket_name (index_num)
		curr_data_set = hdf5_file[chromosome_name+"/"+token_to_search+"/"+bucket_name]

		for curr_chunk in curr_data_set.iter_chunks():
			curr_arr = curr_data_set[curr_chunk]
			if (index_num in curr_arr):
				print ("Found Index: ", end=" ")
				print (index_num)
				return

		print ("Did not Find Index: ", end=" ")
		print (index_num)
		
		hdf5_file.close ()

#*******************************************************************************
#*******************************************************************************

if __name__ == "__main__":

	sq = SeqBuilder ('smallsamplefile', '(?P<CHG>C[ATC]G)|CG|CA|TG', "mytestfile.hdf5")
	sq.process_tokens ()

	sq.check_token ('NC_030850.1', "CAG", 7908)
	sq.check_token ('NC_030850.1', "CHG", 7908)

	sq.check_token ('NC_030850.1', "CAG", 7909)
	sq.check_token ('NC_030850.1', "CHG", 7909)
