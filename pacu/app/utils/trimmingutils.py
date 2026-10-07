import os
from pathlib import Path
from typing import Dict

from pacu import logger


def trimmomatic_dir_adapters() -> Path:
    """
    Retrieves the directory containing the trimmomatic adapters.
    :return: Directory containing trimmomatic adapters
    """
    # Explicitly set adapter directory
    if os.environ.get('TRIMMOMATIC_ADAPTER_DIR') is not None:
        return Path(os.environ['TRIMMOMATIC_ADAPTER_DIR'])

    # CONDA (or Pixi) installation
    if os.environ.get('CONDA_PREFIX') is not None:
        return Path(os.environ['CONDA_PREFIX'], 'share', 'trimmomatic', 'adapters')
    logger.error('TRIMMOMATIC_ADAPTER_DIR environment variable not set')
    raise RuntimeError('TRIMMOMATIC_ADAPTER_DIR environment variable not set')


def trimmomatic_collect_output(dir_out: Path) -> Dict[str, Path]:
    """
    Collects the trimmomatic output.
    :param dir_out: Output directory
    :return: Output dictionary
    """
    try:
        fq_dict: dict[str, Path | None] = {
            '1P': next(dir_out.glob('*_1P.fastq.gz')),
            '2P': next(dir_out.glob('*_2P.fastq.gz'))
        }
    except StopIteration:
        raise FileNotFoundError('No paired end reads survived trimming')

    # Collect output (single-end reads)
    for key in ('1U', '2U'):
        try:
            path_fq = next(dir_out.glob(f'*_{key}.fastq.gz'))
            fq_dict[key] = path_fq
        except StopIteration:
            logger.warning(f"Trimmed unpaired reads ('{key}') not found")
            fq_dict[key] = None
    return fq_dict
