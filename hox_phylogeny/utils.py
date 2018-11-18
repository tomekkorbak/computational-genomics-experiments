from tqdm import tqdm

from Bio import SeqIO, SearchIO
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Phylo import draw_ascii, draw


def prepare_query(content):
	with open('query.txt', 'w') as query_file:
		query_file.write(str(content))


def make_score_matrix(records):
	record_ids = [record.id for record in records]
	matrix = DistanceMatrix(names=record_ids)
	for i, sequence_a in enumerate(tqdm(records)):
		prepare_query(sequence_a.seq)
		blastp_cline = NcbiblastpCommandline(query='query.txt', db='db', outfmt=5, out='./result.txt')
		stdout, stderr = blastp_cline()
		results = SearchIO.read('./result.txt', format='blast-xml')
		for hit in results:
			highest_scoring_pair = max(list(hit), key=lambda hit: hit.bitscore)
			score = highest_scoring_pair.bitscore
			length = len(list(highest_scoring_pair.fragments))
			try:
				j = record_ids.index(highest_scoring_pair.hit_id)
				matrix[i, j] = score/length
			except:
				pass
	return matrix


def make_tree(distance_matrix):
	constructor = DistanceTreeConstructor()
	njtree = constructor.nj(distance_matrix)
	njtree.root_at_midpoint()
	return njtree


def draw_tree(tree, duplications):

	def branch_labels(node):
		if node in duplications:
			return '*'
		else:
			return ''
	
	def label_colors(label):
		mapping = {
			'Mouse': 'red',
			'Human': 'green',
			'Beetle': 'blue',
			'Fruitfly': 'purple'
		}
		for species, color in mapping.items():
			if label.startswith(species):
				return color
		return 'black'

	from pylab import rcParams
	rcParams['figure.figsize'] = 15, 40
	draw(tree, branch_labels=branch_labels, label_colors=label_colors)
