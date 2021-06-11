import tarfile
import requests
import argparse
import zipfile
import tqdm
import tempfile
import re
import shutil
import os
import glob
from gff_longest_transcript import find_longest_transcript
import gzip


def download_file(url, dest=None, description=None, overwrite=False, compressed=False):
    if description is None:
        description = 'Downloading file "{}" ==> "{}"...'.format(url, dest)

    if not overwrite and dest and os.path.isfile(dest):
        print(description)
        return dest

    if re.match("^(http|ftp)", url):
        response = requests.get(url, stream=True)
        total_size_in_bytes= int(response.headers.get('content-length', 0))
        block_size = 1024 #1 Kibibyte

        progress_bar = tqdm.tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
        with tempfile.NamedTemporaryFile(mode="wb", delete=False) as file:
            path = file.name
            progress_bar.set_description(description)
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)
            progress_bar.close()
        if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
            print("ERROR, something went wrong")
    else: # local file
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.close()
        shutil.copy2(url, tmp.name)
        path = tmp.name

    if compressed:
        with gzip.open(path, 'rb') as f_in:
            with tempfile.NamedTemporaryFile(mode="wb", delete=False) as f_out:
                shutil.copyfileobj(f_in, f_out)
                path = f_out.name

    if dest:
        if not os.path.exists(os.path.dirname(dest)):
            os.makedirs(os.path.dirname(dest))

        shutil.copy2(path, dest)
        path = dest

    return path



def download_raw_genome(url, dest, overwrite=False):
    if not overwrite and dest and os.path.isfile(dest):
        print('Downloading raw genome "{}" ==> "{}". Already exists, skipping...'.format(url, dest))
        return

    if not os.path.exists(os.path.dirname(dest)):
        os.makedirs(os.path.dirname(dest))

    path = download_file(url, description='Downloading raw genome "{}"\n'.format(url))

    dest_chromFa = "{}_{}".format(dest, "chromFa")
    with tarfile.open(path, mode="r:gz") as tf:
        progress_bar = tqdm.tqdm(tf.getmembers(), desc="Extracting raw genome archive contents into '{}'\n".format(dest_chromFa))
        for member in progress_bar:
            try:
                tf.extract(member, dest_chromFa)
            except tarfile.TarError as e:
                pass

    print("Joining chromosomes into single fasta file '{}'\n".format(dest))
    with open(dest, 'w') as outfile:
        # Iterate through list
        for names in glob.glob(os.path.join(dest_chromFa, "*")):
            with open(names) as infile:
                outfile.write(infile.read())
            outfile.write("\n")

    shutil.rmtree(dest_chromFa)


def download_bowtie2_index(url, dest, overwrite=False):
    if not overwrite and dest and len(glob.glob(os.path.join(dest, "*.bt2"))) > 0:
        print('Downloading bowtie2 index "{}" ==> "{}". Already exists, skipping...'.format(url, dest))
        return

    path = download_file(url, description="Downloading bowtie index '{}'\n".format(url))

    with zipfile.ZipFile(path, 'r') as zf:
        for member in tqdm.tqdm(zf.infolist(), desc="Extracting bowtie2 index into '{}'\n".format(dest)):
            try:
                zf.extract(member, dest)
            except zipfile.error as e:
                pass

def download_genome(genome, path):
    download_file("http://hgdownload.cse.ucsc.edu/goldenpath/{genome}/database/chromInfo.txt.gz".format(genome=genome), os.path.join(path, "{genome}/annotation/ChromInfo.txt".format(genome=genome)), compressed=True)
    download_file("http://hgdownload.cse.ucsc.edu/goldenpath/{genome}/database/cytoBand.txt.gz".format(genome=genome), os.path.join(path, "{genome}/annotation/cytoBand.txt".format(genome=genome)), compressed=True)
    download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/chromFa.tar.gz".format(genome=genome), os.path.join(path, "{genome}/{genome}.fa".format(genome=genome)))

    download_file("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.chrom.sizes".format(genome=genome), os.path.join(path, "{genome}/annotation/{genome}.chrom.sizes".format(genome=genome)))
    download_bowtie2_index("https://genome-idx.s3.amazonaws.com/bt/{genome}.zip".format(genome=genome), path)
    #
    download_file("http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/genes/{genome}.refGene.gtf.gz".format(genome=genome), os.path.join(path, "{genome}.refGene.gtf.gz".format(genome=genome)))
    print("Creating annotation file from {gtf}...".format(gtf=os.path.join(path, "{genome}.refGene.gtf.gz".format(genome=genome))))
    find_longest_transcript(os.path.join(path, "{genome}.refGene.gtf.gz".format(genome=genome)), os.path.join(path, "{genome}/annotation/refGene.bed".format(genome=genome)), clip_start=50, clip_strand_specific=True)


def download_dependencies(path):
    # # Download libraries used in the pipeline
    download_file("https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip", os.path.join(path, "bowtie2-2.2.9-linux-x86_64.zip"))
    download_file("https://master.dl.sourceforge.net/project/samtools/samtools/1.6/samtools-1.6.tar.bz2", os.path.join(path, "samtools-1.6.tar.bz2"))
    download_file("https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz", os.path.join(path, "bedtools-2.29.1.tar.gz"))
    download_file("https://anaconda.org/bioconda/gfold/1.1.4/download/linux-64/gfold-1.1.4-gsl2.2_2.tar.bz2", os.path.join(path, "gfold-1.1.4-gsl2.2_2.tar.bz2"))
    download_file("http://www.artfiles.org/gnu.org/gsl/gsl-2.2.1.tar.gz", os.path.join(path, "gsl-2.2.1.tar.gz"))
    download_file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig", os.path.join(path, "wigToBigWig"))
    download_file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig", os.path.join(path, "bigWigToWig"))
    download_file("https://github.com/bxlab/bx-python/archive/v0.8.9.tar.gz", os.path.join(path, "bx-python-0.8.9.tar.gz"))
    download_file("https://netcologne.dl.sourceforge.net/project/rseqc/RSeQC-2.6.4.tar.gz", os.path.join(path, "RSeQC-2.6.4.tar.gz"))
    download_file("https://github.com/numpy/numpy/releases/download/v1.16.6/numpy-1.16.6.tar.gz", os.path.join(path, "numpy-1.16.6.tar.gz"))
    download_file("https://github.com/cython/cython/archive/0.29.19.tar.gz", os.path.join(path, "cython-0.29.19.tar.gz"))
    download_file("https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/python-nose/nose-1.0.0%20(1).tar.gz", os.path.join(path, "nose-1.0.0.tar.gz"))
    download_file("https://github.com/pysam-developers/pysam/archive/v0.16.0.tar.gz", os.path.join(path, "pysam-0.16.0.tar.gz"))
    download_file("https://github.com/benjaminp/six/archive/1.15.0.tar.gz", os.path.join(path, "six-1.15.0.tar.gz"))
    download_file("https://github.com/psf/requests/archive/v2.24.0.tar.gz", os.path.join(path, "requests-2.24.0.tar.gz"))
    download_file("https://github.com/tqdm/tqdm/archive/v4.49.0.tar.gz", os.path.join(path, "tqdm-4.49.0.tar.gz"))
    download_file("https://github.com/daler/pybedtools/archive/v0.8.0.tar.gz", os.path.join(path, "pydevtools-0.8.0.tar.gz"))
    download_file("https://github.com/pandas-dev/pandas/archive/v0.24.2.tar.gz", os.path.join(path, "pandas-0.24.2.tar.gz"))
    download_file("https://github.com/dateutil/dateutil/releases/download/2.8.1/python-dateutil-2.8.1.tar.gz", os.path.join(path, "dateutil-2.8.1.tar.gz"))
    download_file("https://github.com/kjd/idna/archive/v2.10.tar.gz", os.path.join(path, "idna-2.10.tar.gz"))
    download_file("https://github.com/urllib3/urllib3/archive/1.25.10.tar.gz", os.path.join(path, "urllib3-1.25.10"))
    download_file("https://github.com/chardet/chardet/archive/3.0.4.tar.gz", os.path.join(path, "chardet-3.0.4.tar.gz"))
    download_file("https://github.com/certifi/python-certifi/archive/2020.06.20.tar.gz", os.path.join(path, "python-certifi-2020.06.20.tar.gz"))



    #download_file("https://github.com/stub42/pytz/archive/release_2020.1.tar.gz", os.path.join(path, "pytz-2020.1.tar.gz"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download groseq dependencies')
    parser.add_argument('data', choices=['mm9', 'mm10', 'hg19', 'hg38', 'dependencies'], help="""Should be one of the following: \n
    mm9 - mm9 model files\n  
    mm10 - mm10 model files\n
    hg19 - hg19 model files""")
    parser.add_argument('path', nargs="?", default=".", help="Path where the files will be downloaded (default: current directory)")

    args = parser.parse_args()

    if not os.path.exists(args.path):
        os.makedirs(args.path)

    if args.data == "dependencies":
        download_dependencies(args.path)
    else:
        download_genome(args.data, args.path)