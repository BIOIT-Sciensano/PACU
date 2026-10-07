import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from pacu.app.utils.loggingutils import initialize_logging
from pacu.map_to_ref import MapToRef


class TestMapToRef(unittest.TestCase):
    """
    Tests the map to ref script.
    """

    def test_map_to_ref_illumina(self) -> None:
        """
        Tests the map to ref script with Illumina data.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--read-type', 'illumina',
                '--fastq-illumina',
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_1P.fastq.gz')),
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_2P.fastq.gz')),
                '--output', str(path_out),
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_map_to_ref_with_fasta_name(self) -> None:
        """
        Tests the map to ref script with the fasta name parameter set.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--ref-fasta-name', 'my_reference.fasta',
                '--read-type', 'illumina',
                '--fastq-illumina',
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_1P.fastq.gz')),
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_2P.fastq.gz')),
                '--output', str(path_out),
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_map_to_ref_illumina_with_trim(self) -> None:
        """
        Tests the map to ref script with Illumina data and trimming enabled.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--read-type', 'illumina',
                '--fastq-illumina',
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_1P.fastq.gz')),
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_2P.fastq.gz')),
                '--output', str(path_out),
                '--trim',
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_map_to_ref_ont(self) -> None:
        """
        Tests the map to ref script with ONT data.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--read-type', 'ont',
                '--fastq-ont',
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151-ont.fastq.gz')),
                '--output', str(path_out),
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_map_to_ref_ont_with_trim(self) -> None:
        """
        Tests the map to ref script with ONT data and trimming enabled.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--read-type', 'ont',
                '--fastq-ont',
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151-ont.fastq.gz')),
                '--trim',
                '--output', str(path_out),
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_map_to_ref_rerun_same_working_dir(self) -> None:
        """
        Tests the map to ref script when it is executed twice in the same working directory (incl. Galaxy names).
        :return: None
        """
        dir_testdata = files('pacu').joinpath('resources/testdata')
        args_by_read_type = {
            'ont': [
                '--fastq-ont', str(dir_testdata.joinpath('fastq/TIAC1151-ont.fastq.gz')),
                '--fastq-ont-name', 'galaxy_ont.fastq.gz'
            ],
            'illumina': [
                '--fastq-illumina',
                str(dir_testdata.joinpath('fastq/TIAC1151_1P.fastq.gz')),
                str(dir_testdata.joinpath('fastq/TIAC1151_2P.fastq.gz')),
                '--fastq-illumina-names', 'galaxy_ilmn_1P.fastq.gz', 'galaxy_ilmn_2P.fastq.gz'
            ]
        }
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            for read_type, args_fastq in args_by_read_type.items():
                path_out = Path(dir_, f'mapped_reads_{read_type}.bam')
                for _ in range(2):
                    map_to_ref = MapToRef([
                        '--ref-fasta', str(dir_testdata.joinpath('NC_002695.2-subset.fasta')),
                        '--read-type', read_type,
                        *args_fastq,
                        '--output', str(path_out),
                        '--dir-working', dir_,
                        '--threads', '4'
                    ])
                    map_to_ref.run()

                # Verify the output file
                self.assertTrue(path_out.exists())
                self.assertGreater(path_out.stat().st_size, 0)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
