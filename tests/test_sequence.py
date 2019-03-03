import unittest
import os

import sys
sys.path.append(os.getenv('PWD'))
import sequencer

class TestFilterSequence(unittest.TestCase):

	def setUp(self):
		pass

	def test_uppercase_input(self):
		nucleotides = ['A', 'G', 'T', 'C', 'G']
		scores = ['10', '19', '50', '40', '30']
		filtered_sequence = sequencer.FilterSequence(nucleotides, scores, 20)
		self.assertEqual(filtered_sequence, ['-', '-', 'T', 'C', 'G'])

	def test_lowercase_input(self):
		nucleotides = ['a', 'g', 't', 'c', 'g']
		scores = ['10', '19', '50', '40', '30']
		filtered_sequence = sequencer.FilterSequence(nucleotides, scores, 20)
		self.assertEqual(filtered_sequence, ['-', '-', 't', 'c', 'g'])

	def test_mixcase_input(self):
		nucleotides = ['a', 'G', 'T', 'c', 'g']
		scores = ['20', '19', '50', '40', '5']
		filtered_sequence = sequencer.FilterSequence(nucleotides, scores, 20)
		self.assertEqual(filtered_sequence, ['a', '-', 'T', 'c', '-'])

	def test_bad_input(self):
		nucleotides = ['y', 'g', 'E', 'R', 'R']
		scores = ['g', 'g', 'h', '&', 'x']
		self.assertRaises(ValueError, sequencer.FilterSequence, nucleotides, scores, 20)

	def test_single_long_string(self):
		nucleotides = 'AGCAA'
		scores = "20 19 5 30 5"
		filtered_sequence = sequencer.FilterSequence(list(nucleotides), scores.split(' '), 20)
		self.assertEqual(filtered_sequence, ['A', '-', '-', 'A', '-'])

	def test_unequal_datasets(self):
		nucleotides = 'AGCA'
		scores = "20 30 40"
		with self.assertRaises(ValueError):
			sequencer.FilterSequence(list(nucleotides), scores.split(' '), 20)

		#The following code won't work, left in here as a historical reminder. This is because
		#unittest's assertRaises needs a wrapper to be able to access the raised exception.
		#self.assertRaises(ValueError, sequencer.FilterSequence, nucleotides, scores, 20)


class TestReverseComplement(unittest.TestCase):
	def setUp(self):
		pass

	def test_uppercase_input(self):
		sequence = "AGTCAT"
		self.assertEqual(sequencer.ReverseComplement(sequence), "ATGACT")

	def test_mixedcase_input(self):
		sequence = "agtCAt"
		self.assertEqual(sequencer.ReverseComplement(sequence), "aTGact")

	def test_lowercase_input(self):
		sequence = "agtcat"
		self.assertEqual(sequencer.ReverseComplement(sequence), "atgact")

	def test_gaps_input(self):
		sequence = "agt-C"
		self.assertEqual(sequencer.ReverseComplement(sequence), "G-act")


if __name__ == '__main__':
    unittest.main()