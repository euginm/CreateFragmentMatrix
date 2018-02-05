class FragMatrixCreator:
    def __init__(self, bam_path, vcf_path, region, output_folder, output_prefix):
        from pathlib import Path
        import pysam
        import subprocess

        self.region_name = region
        print('*** Processing region ' + self.region_name + '... ***', sep='')

        with pysam.AlignmentFile(bam_path) as bam:
            try:
                bam.check_index()
            except ValueError:
                print('*** ' + bam_path + '.bai not found, indexing bam file... ***')
                subprocess.check_call(['samtools', 'index', bam_path])

        self.vcf_path = Path(vcf_path).resolve()
        self.vcf_index_path = Path(vcf_path + '.index')
        if not self.vcf_index_path.exists():
            print('*** ' + str(self.vcf_index_path) + ' doesn\'t exist, indexing vcf... ***')
            self.index_vcf_file()

        self.vcf_position = self.get_vcf_position()

        output_folder = Path(output_folder).resolve()
        try:
            output_folder.mkdir(parents=True)
        except FileExistsError:
            pass
        if output_prefix:
            output_prefix += '_'
        self.output_matrix_path = output_folder / (output_prefix + self.region_name + '.frags')
        self.output_genotypes_path = output_folder / (output_prefix + self.region_name + '.genotypes')

        self.polymorphisms = {}  # {position: [ref, alt1, alt2, ...]}

        bam = pysam.AlignmentFile(bam_path)
        self.alignments = bam.fetch(self.region_name)

    def index_vcf_file(self):
        """
        Index a vcf file, write it to vcf_file.index
            region name         position in vcf file
        """
        
        with self.vcf_path.open() as vcf_file, self.vcf_index_path.open('w') as index_file:
            #  skip headers and store pointer position
            while True:
                line = vcf_file.readline()
                if line[:2] != '##':
                    vcf_position = vcf_file.tell()
                    line = vcf_file.readline()
                    break

            while line:
                region_name = line.split('\t')[0]
                index_file.write(region_name + '\t' + str(vcf_position) + '\n')

                while line and region_name == line.split('\t')[0]:
                    vcf_position = vcf_file.tell()
                    line = vcf_file.readline()
        return

    def get_vcf_position(self):
        """
        :return: INT, position of genomic region in vcf file
        """
        import sys

        with self.vcf_index_path.open() as vcf_index:
            for line in vcf_index:
                region_name = line.split()[0]
                position_in_vcf = line.split()[1]
                if region_name == self.region_name:
                    return int(position_in_vcf)
        sys.exit('ERROR ' + self.region_name + ' is not found in ' + str(self.vcf_index_path) + ' ERROR')

    def get_region_polymorphisms(self):
        """
        Fill the region_polymorphisms dictionary
        """
        import sys

        print('*** Loading polymorphic sites from ', self.vcf_path, '... ***', sep='')
        with self.vcf_path.open() as vcf_file:
            vcf_file.seek(self.vcf_position)
            for line in vcf_file:
                line = line.split()
                if line[0] != self.region_name:
                    print('*** Loading polymorphic sites done! ***')
                    return
                self.polymorphisms[int(line[1]) - 1] = [line[3]] + line[4].split(',')  # -1 for compatibility with pysam

        if len(self.polymorphisms) < 2:
            sys.exit('ERROR ' + self.region_name + ' has less than 2 polymorphic sites! ERROR')
        print('*** Loading polymorphic sites done! ***')
        return

    def write_fragment_matrix(self):
        print('*** Converting alignments to fragments... ***')

        # Enumerate polymorphism positions and map them to their index
        sorted_polymorphism_positions = sorted(self.polymorphisms.keys())
        polymorphism_to_index = {pos: str(index) for index, pos in enumerate(sorted_polymorphism_positions, start=1)}

        fragments = []  # [read name, pos, alleles, qualities, mapq]
        # store fragment names which have only one allele, write them to matrix
        # if mate has at least one informative allele
        fragments_purgatorium = {}

        for alignment in self.alignments:
            if alignment.is_unmapped:
                continue

            # create the read map with {reference_pos: abs_read_pos}
            alignment_map = {ref_pos: read_pos for read_pos, ref_pos in alignment.get_aligned_pairs(matches_only=True)}

            # store alignment and polymorphism overlapping positions
            fragment_positions = list(range(alignment.reference_start, alignment.reference_end) &
                                      self.polymorphisms.keys())

            # fill the fragment dictionary with alleles: {polymorphic_site: allele number}
            query_sequence = alignment.query_sequence
            query_qualities = alignment.query_qualities
            alignment_name = alignment.query_name
            fragment = {}
            qualities = {}

            for pos in fragment_positions:
                allele = ''
                polymorphism_length = len(self.polymorphisms[pos][0])
                qual = 0
                if polymorphism_length == 1:  # in case of SNP
                    if pos in alignment_map:
                        allele = query_sequence[alignment_map[pos]]
                        qual = query_qualities[alignment_map[pos]]
                else:  # in case of complex polymorphism
                    if pos in alignment_map and pos + polymorphism_length - 1 in alignment_map:
                        if pos + polymorphism_length not in alignment_map:
                            allele = query_sequence[alignment_map[pos]:]
                            quals = query_qualities[alignment_map[pos]:]
                        else:
                            allele = query_sequence[alignment_map[pos]:alignment_map[pos + polymorphism_length]]
                            quals = query_qualities[alignment_map[pos]:alignment_map[pos + polymorphism_length]]
                        qual = round(sum(quals) / len(quals))

                if allele not in self.polymorphisms[pos]:
                    continue
                fragment[pos] = str(self.polymorphisms[pos].index(allele))
                qualities[pos] = chr(qual + 33)  # Phred+33

            if not fragment:  # no alleles
                continue

            # add positions without alleles between fragment start and end
            fragment_positions = sorted(list(range(min(fragment.keys()), max(fragment.keys()) + 1) &
                                             self.polymorphisms.keys()))
            fragment_string = ''
            qualities_string = ''
            for pos in fragment_positions:
                fragment_string += fragment[pos] if pos in fragment else '-'
                qualities_string += qualities[pos] if pos in qualities else '-'

            if len(fragment) == 1:  # only one allele decisions
                if alignment_name in fragments_purgatorium:  # a mate with one allele is already present
                    # check if allele positions are not overlapping
                    if fragments_purgatorium[alignment_name] != polymorphism_to_index[min(fragment.keys())]:
                        del fragments_purgatorium[alignment_name]
                    else:
                        continue
                elif alignment.is_paired and alignment.mate_is_unmapped:  # one allele and the mate is unmapped
                    continue
                elif alignment.next_reference_name != self.region_name:  # the mate is mapped to another region
                    continue
                else:  # else add to waiting list
                    fragments_purgatorium[alignment_name] = polymorphism_to_index[min(fragment.keys())]
            elif alignment_name in fragments_purgatorium:  # a mate with one allele present, delete it from waiting list
                del fragments_purgatorium[alignment_name]

            fragments.append([alignment_name, polymorphism_to_index[min(fragment.keys())], fragment_string,
                              qualities_string, str(alignment.mapping_quality)])

        if not fragments:
            print('*** No informative fragments! ***')
            return

        print('*** Writing ', self.output_matrix_path, '... ***', sep='')
        with self.output_matrix_path.open('w') as output_matrix:
            output_matrix.write('>' + self.region_name + '\n')
            output_matrix.write('\n'.join(['\t'.join(fragment_representation) for fragment_representation in fragments
                                           if fragment_representation[0] not in fragments_purgatorium]))
        print('*** Writing ', self.output_matrix_path, ' done! ***', sep='')
        return

    def write_genotypes(self):
        output_genotypes = self.output_genotypes_path.open('w')
        pass
        return

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Convert bam and vcf into fragment matrix')
    parser.add_argument('bam_path', help='Path to bam file')
    parser.add_argument('vcf_path', help='Path to vcf file')
    parser.add_argument('region', help='Genomic region name')
    parser.add_argument('output_folder', help='Output folder name')
    parser.add_argument('--output_prefix', help='Specify the output file prefix')
    parser.add_argument('--genotypes', action='store_true',
                        help='Set this flag to create the genotypes file along with fragments matrix')
    args = vars(parser.parse_args())
    if args.get('output_prefix') is None:
        args['output_prefix'] = ''
    fragMatrixCreator = FragMatrixCreator(args.get('bam_path'), args.get('vcf_path'), args.get('region'),
                                          args.get('output_folder'), args.get('output_prefix'))
    fragMatrixCreator.get_region_polymorphisms()
    fragMatrixCreator.write_fragment_matrix()
    if args.get('genotypes'):
        fragMatrixCreator.write_genotypes()
