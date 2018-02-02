class FragMatrixCreator:
    def __init__(self, bam_path, vcf_path, region, output_folder, output_prefix='', genotypes=False):
        import os.path
        import pysam
        import subprocess

        self._bam_path = os.path.abspath(bam_path)
        with pysam.AlignmentFile(bam_path) as bam:
            try:
                bam.check_index()
            except ValueError:
                print '*** ' + bam_path + '.bai not found, indexing bam file... ***'
                subprocess.check_call(['samtools', 'index', bam_path])
                
        self._vcf_path = os.path.abspath(vcf_path)
        self._vcf_index_path = vcf_path + '.index'
        if not os.access(self._vcf_index_path, os.F_OK):
            print '*** ' + self._vcf_index_path + ' doesn\'t exist, indexing vcf... ***'
            self.index_vcf_file()

        self._region_name = region
        self._vcf_position = self.get_vcf_position()
        
        self._output_folder = os.path.abspath(output_folder)
        self._output_matrix_path = os.path.join(self._output_folder, output_prefix + self._region_name + '.frags')
        
        self._create_genotypes = genotypes
        if self._create_genotypes:
            self._output_genotypes_path = os.path.join(self._output_folder, output_prefix + self._region_name + '.genotypes')

        self._region_polymorphisms = {}  # {position: [ref, alt1, alt2, ...]}

    def index_vcf_file(self):
        """
        Index a vcf file, write it to vcf_file.index
            region name         position in vcf file
        """
        
        with open(self._vcf_path, 'r') as vcf_file, open(self._vcf_index_path, 'w') as index_file:
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
        with open(self._vcf_index_path) as vcf_index:
            for line in vcf_index:
                region_name = line.split()[0]
                position_in_vcf = line.split()[1]
                if region_name == self._region_name:
                    return int(position_in_vcf)
        sys.exit('ERROR ' + self._region_name + ' is not found in ' + self._vcf_index_path + ' ERROR')

    def get_region_polymorphisms(self):
        """
        Fill the region_polymorphisms dictionary
        """
        with open(self._vcf_path) as vcf_file:
            vcf_file.seek(self._vcf_position)
            for line in vcf_file:
                line = line.split()
                if line[0] != self._region_name:
                    return

                self._region_polymorphisms[int(line[1]) - 1] = [line[3]] + line[4].split(',')  # -1 for compatibility with pysam
        return

    def write_fragment_matrix(self):
        pass

    def write_genotypes(self):
        pass

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
