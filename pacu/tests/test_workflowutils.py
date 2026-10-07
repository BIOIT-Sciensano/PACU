import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

import pandas as pd

from pacu.app.utils import workflowutils


class TestWorkflowUtils(unittest.TestCase):
    """
    Tests for the workflow utils.
    """

    def test_plot_newick_phylogeny(self) -> None:
        """
        Tests the plot_newick_phylogeny function.
        """
        path_nwk = Path(str(files('pacu').joinpath('resources/testdata/phylogeny.nwk')))
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            # Create visualization
            path_out = Path(dir_, 'tree.png')
            workflowutils.plot_newick_phylogeny(path_nwk, path_out)

            # Verify output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_determine_name_from_fq_ont(self) -> None:
        """
        Tests the 'determine_name_from_fq' function with ONT input.
        :return: None
        """
        self.assertEqual(workflowutils.determine_name_from_fq(fq_ont=Path('/path/to/reads.fastq')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(fq_ont=Path('/path/to/reads.fastq.gz')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(fq_ont=Path('/path/to/reads.fq')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(fq_ont=Path('/path/to/reads.fq.gz')), 'reads')

    def test_determine_name_from_fq_illumina(self) -> None:
        """
        Tests the 'determine_name_from_fq' function with Illumina input.
        :return: None
        """
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/my-reads_1P.fastq')), 'my-reads')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/input_file_1.fastq.gz')), 'input_file')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_1.fq')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_1P.fq.gz')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_S22_L001_R1_001.fastq.gz')), 'reads_S22')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_S22_L001_R1_001.fq.gz')), 'reads_S22')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/unknown_pattern.fastq')), 'unknown_pattern')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/unknown_pattern.fastq.gz')), 'unknown_pattern')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_R1.fq')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_R1.fastq.gz')), 'reads')

    def test_parse_bed(self) -> None:
        """
        Tests the 'parse_bed' function with empty files, header lines and additional columns.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_empty = Path(dir_, 'empty.bed')
            path_empty.touch()
            self.assertEqual(len(workflowutils.parse_bed(path_empty)), 0)
            self.assertEqual(workflowutils.count_regions(path_empty), 0)
            self.assertEqual(workflowutils.count_covered_positions(path_empty), 0)

            path_bed = Path(dir_, 'regions.bed')
            path_bed.write_text('track name=phages\n# comment\n\nchr1\t10\t20\tphage_1\t0\t+\n2\t5\t8\n')
            data_bed = workflowutils.parse_bed(path_bed)
            self.assertEqual(data_bed.to_dict('records'), [
                {'chr': 'chr1', 'start': 10, 'end': 20},
                {'chr': '2', 'start': 5, 'end': 8}
            ])
            self.assertEqual(workflowutils.count_regions(path_bed), 2)
            self.assertEqual(workflowutils.count_regions(str(path_bed)), 2)
            self.assertEqual(workflowutils.count_covered_positions(path_bed), 13)

    def test_is_new_region(self) -> None:
        """
        Tests the 'is_new_region' function, including consecutive positions on different contigs.
        :return: None
        """
        data = pd.DataFrame({'chr': ['c1', 'c1', 'c1', 'c2', 'c2'], 'pos': [1, 2, 5, 6, 7]})
        data['shift_pos'] = data['pos'].shift()
        data['shift_chr'] = data['chr'].shift()
        self.assertEqual(
            [workflowutils.is_new_region(row) for _, row in data.iterrows()], [False, False, True, True, False])

    def test_calculate_overlaps(self) -> None:
        """
        Tests the 'calculate_overlaps' function with a multi-contig reference.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_phages = Path(dir_, 'phages.bed')
            path_phages.write_text('c1\t0\t3\n')
            path_gubbins = Path(dir_, 'gubbins.bed')
            path_gubbins.write_text('c2\t3\t5\n')
            path_depth = Path(dir_, 'depth.bed')
            path_depth.write_text('c1\t2\t4\n')
            data_overlap = workflowutils.calculate_overlaps(
                {'c1': 10, 'c2': 5}, path_phages, path_gubbins, path_depth)

        self.assertEqual(data_overlap.to_dict('records'), [
            {'chr': 'c1', 'pos': 1, 'phages': True, 'gubbins': False, 'depth': False},
            {'chr': 'c1', 'pos': 2, 'phages': True, 'gubbins': False, 'depth': False},
            {'chr': 'c1', 'pos': 3, 'phages': True, 'gubbins': False, 'depth': True},
            {'chr': 'c1', 'pos': 4, 'phages': False, 'gubbins': False, 'depth': True},
            {'chr': 'c2', 'pos': 4, 'phages': False, 'gubbins': True, 'depth': False},
            {'chr': 'c2', 'pos': 5, 'phages': False, 'gubbins': True, 'depth': False},
        ])

    def test_sanitize_input_name(self) -> None:
        """
        Tests the 'sanitize_input_name' function.
        :return: None
        """
        self.assertEqual(workflowutils.sanitize_input_name('my sample.bam', 'bam'), 'my_sample.bam')
        self.assertEqual(workflowutils.sanitize_input_name('my sample', 'bam'), 'my_sample.bam')
        self.assertEqual(workflowutils.sanitize_input_name('a/b#c.bam', 'bam'), 'abc.bam')
        self.assertEqual(workflowutils.sanitize_input_name('a/b#c', 'bam'), 'abc.bam')
        self.assertEqual(workflowutils.sanitize_input_name('sample.', 'fasta'), 'sample.fasta')


if __name__ == '__main__':
    unittest.main()
