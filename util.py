import random
from logging import StreamHandler, getLogger as realGetLogger, Formatter
import numpy as np
import pandas as pd
from colorama import Fore, Back, Style


###############################################################################
# Configure color stream handler
# https://gist.github.com/jonaprieto/a61d9cade3ba19487f98
###############################################################################
class ColourStreamHandler(StreamHandler):
    """ A colorized output StreamHandler """

    # Some basic colour scheme defaults
    colours = {
        'DEBUG': Fore.CYAN,
        'INFO': Fore.GREEN,
        'WARN': Fore.YELLOW,
        'WARNING': Fore.YELLOW,
        'ERROR': Fore.RED,
        'CRIT': Back.RED + Fore.WHITE,
        'CRITICAL': Back.RED + Fore.WHITE
    }

    def emit(self, record):
        try:
            message = self.format(record)
            self.stream.write(self.colours[record.levelname] + message +
                              Style.RESET_ALL)
            self.stream.write(getattr(self, 'terminator', '\n'))
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


###############################################################################
# configure logger
###############################################################################
def get_logger(name=None,
               fmt='[%(name)s::%(funcName)s] %(levelname)s %(message)s',
               level='INFO'):
    """ Get and initialize a colourised logging instance if the system supports
    it as defined by the log.has_colour
    :param name: Name of the logger
    :type name: str
    :param fmt: Message format to use
    :type fmt: str
    :return: Logger instance
    :rtype: Logger
    """
    log = realGetLogger(name)
    # Only enable colour if support was loaded properly
    handler = ColourStreamHandler()
    handler.setLevel(level)
    handler.setFormatter(Formatter(fmt))
    log.addHandler(handler)
    log.setLevel(level)
    log.propagate = 0  # Don't bubble up to the root logger
    return log


util_log = get_logger(__name__, level="DEBUG")


###############################################################################
# process simulated data
###############################################################################
def process_ts(variants, model_label, rep_label, mnm_frac, mnm_dist, method,
               out_dir):
    if method == "archie":
        eig_data = generateEigData(
            variants,
            model_label,
            rep_label=rep_label,
            mnm_frac=mnm_frac,
            mnm_dist=mnm_dist)

        util_log.debug("---Generating data without MNMs---")
        util_log.debug("Prefix: %s", eig_data.prefix)
        util_log.debug("Genotype matrix dimensions: %s, %s" % tuple(
            pd.DataFrame(eig_data.geno).shape))
        util_log.debug("SNP table dimensions: %s, %s" % tuple(
            pd.DataFrame(eig_data.snp).shape))
        writeEigData(eig_data, out_dir).dump()

        eig_data_mnms = eig_data.mnms()
        util_log.debug("---Generating data with MNMs---")
        util_log.debug("Prefix: %s", eig_data_mnms.prefix)
        util_log.debug("Genotype matrix dimensions: %s, %s" % tuple(
            eig_data_mnms.geno.shape))
        util_log.debug(
            "SNP table dimensions: %s, %s" % tuple(eig_data_mnms.snp.shape))
        writeEigData(eig_data_mnms, out_dir).dump()

    elif method == "sprime":
        vcf_data = generateVCFData(
            variants,
            model_label,
            rep_label=rep_label,
            mnm_frac=mnm_frac,
            mnm_dist=mnm_dist)

        util_log.debug("---Generating data without MNMs---")
        util_log.debug("Prefix: %s", vcf_data.prefix)
        util_log.debug(
            "VCF dimensions: %s, %s" % tuple(pd.DataFrame(vcf_data.vcf).shape))
        writeVCFData(vcf_data, out_dir).dump()

        vcf_data_mnms = vcf_data.mnms()
        util_log.debug("---Generating data with MNMs---")
        util_log.debug("Prefix: %s", vcf_data_mnms.prefix)
        util_log.debug("VCF dimensions: %s, %s" % tuple(
            pd.DataFrame(vcf_data_mnms.vcf).shape))
        writeVCFData(vcf_data_mnms, out_dir).dump()


###############################################################################
# classes for generating data in VCF format
###############################################################################
class mnmVCFData:
    """
    object defines .vcf matrices
    """

    def __init__(self, vcf, prefix, repeat_rows):
        self.vcf = vcf
        self.prefix = prefix
        self.repeat_rows = repeat_rows


class generateVCFData:
    """
    generate VCF files from tree sequence
    """

    def __init__(self,
                 variants,
                 model_label,
                 rep_label=0,
                 mnm_frac=0,
                 mnm_dist=0):
        self.variants = variants
        self.model_label = model_label
        self.rep_label = rep_label
        if mnm_frac > 0:
            self.sim_mnms = True
        else:
            self.sim_mnms = False
        self.mnm_frac = mnm_frac
        self.mnm_dist = mnm_dist
        self.vcf = self.vcf()
        self.prefix = self.prefix()

    def prefix(self):
        """
        get unique identifier for model based on inputs
        """

        prefix = self.model_label + "_rep" + str(self.rep_label)
        return prefix

    def vcf(self):
        """
        generate data frame of genotypes in VCF format
        """

        # for variant in self.ts.variants():
        #     pos = round(variant.site.position)
        # output current variant to VCF

        # text_file.write("\t".join(snv) + "\n")

        snvs = []
        # d = np.empty([59, 2], dtype=int)
        for variant in self.variants:
            pos = round(variant.site.position)
            snv = [
                "1",
                str(pos),
                "1_%s" % pos, "A", "G", "50", "PASS", "VT=SNP", "GT"
            ]
            snv.extend([("|").join(
                [str(g) for g in variant.genotypes[i:i + 2]])
                        for i in range(0, len(variant.genotypes), 2)])
            snvs.append(snv)


        # d = d.append(pd.DataFrame(snv, index=[0]), ignore_index=True)
        # snv_df = pd.DataFrame(snvs)
        # util_log.debug(tuple(pd.DataFrame(snvs).shape))

        return snvs

    def mnms(self):
        """
        update .vcf output to include simulated MNMs
        """

        # snv_df = pd.DataFrame()
        snv_df = []
        #         geno_mat = self.geno()
        repeat_rows = []
        for snv in self.vcf:
            snv_df.append(snv)

            random.seed(int(snv[1]))
            if random.random() < self.mnm_frac:
                # make sure a mnm is at least 10bp from its counterpart.
                # this is to avoid the internal filter of Sprime.
                mnm_snv = snv
                dist = random.randint(10, self.mnm_dist)
                mnm_cand = int(snv[1]) + dist
                mnm_cand_r = str(round(mnm_cand))
                mnm_snv[1] = str(mnm_cand_r)
                mnm_snv[2] = "1_%s" % mnm_cand_r
                snv_df.append(mnm_snv)
                repeat_rows.append(2)
            else:
                repeat_rows.append(1)

        # util_log.debug("\t".join(snv_df[i]) for i in range(0, 4))
        mnm_prefix = self.prefix + "_mnm" + str(self.mnm_dist) + "-" + str(
            self.mnm_frac)

        return mnmVCFData(
            vcf=snv_df, prefix=mnm_prefix, repeat_rows=repeat_rows)


class writeVCFData:
    """
    functions for writing .vcf files
    """

    def __init__(self, vcf_data, output_dir="./"):
        self.vcf_data = vcf_data
        self.prefix = vcf_data.prefix
        self.output_dir = output_dir
        self.out_file = self.output_dir + self.prefix + ".vcf"

    def dump(self):
        """
        run all 3 output functions: write_snp, write_geno, and write_ind
        """

        self.write_vcf()
        # self.write_geno()
        # self.write_ind()

    def write_vcf(self):
        """
        write .vcf files
        """

        header = [
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT"
        ]

        l_IDs = []
        for popID in ["AFR", "EUR", "EA"]:
            #        dict_IDs[popID] = []
            for indID in range(50):
                #            dict_IDs[popID].append("%s_ind%s" % (popID, indID))
                l_IDs.append("%s_ind%s" % (popID, indID))

        header.extend(l_IDs)
        vcf_out = pd.DataFrame(self.vcf_data.vcf)
        colnames = dict(zip(range(0, len(header)), header))
        vcf_out.rename(columns=colnames, inplace=True)

        vcf_out_file = self.output_dir + self.prefix + ".vcf"
        util_log.debug("writing VCF to %s", vcf_out_file)
        vcf_out.to_csv(vcf_out_file, index=False, sep="\t")


###############################################################################
# classes for generating data in EIGENSTRAT format
###############################################################################


class mnmEigData:
    """
    object defines .geno and .snp matrices
    """

    def __init__(self, geno, snp, prefix, repeat_rows):
        self.geno = geno
        self.snp = snp
        self.prefix = prefix
        self.repeat_rows = repeat_rows


class generateEigData:
    """
    class describes the .snp, .geno, and .ind EIGENSTRAT formats
    """

    def __init__(self,
                 variants,
                 model_label,
                 rep_label=0,
                 mnm_frac=0,
                 mnm_dist=0):
        # self.ts = ts
        self.variants = variants
        self.model_label = model_label
        self.rep_label = rep_label
        if mnm_frac > 0:
            self.sim_mnms = True
        else:
            self.sim_mnms = False
        self.mnm_frac = mnm_frac
        self.mnm_dist = mnm_dist
        self.geno = self.geno()
        self.snp = self.snp()
        self.prefix = self.prefix()

    def prefix(self):
        """
        get unique identifier for model based on inputs
        """

        prefix = self.model_label + "_rep" + str(self.rep_label)
        #         if self.sim_mnms:
        #             prefix = prefix + "_mnm"+str(self.mnm_dist)+"-"+str(self.mnm_frac)

        return prefix

    def geno(self):
        """
        generate .geno file from tree sequence
        """

        # util_log.debug(len(self.variants[0].genotypes))
        geno = np.zeros((0, len(self.variants[0].genotypes)), dtype=np.int8)
        # geno = np.array([])
        # idx = 0
        for variant in self.variants:
            # util_log.debug(variant.genotypes)
            # util_log.debug(geno)
            # util_log.debug(geno[idx])
            geno = np.append(geno, np.asarray([variant.genotypes]), axis=0)
            # idx += 1

        # geno = np.delete(geno, (0), axis=0)  # remove dummy 1st row
        return geno

    def snp(self):
        """
        generate .snp file from tree sequence
        """

        d = pd.DataFrame(columns=['ID', 'CHR', 'POS1', 'POS', 'REF', 'ALT'])
        # for variant in list(self.ts.variants()):
        for variant in self.variants:
            d = d.append(
                pd.DataFrame({
                    'ID': "1:" + str(round(variant.site.position)),
                    'CHR': "1",
                    'POS1': str(variant.site.position / 10e6),
                    'POS': str(round(variant.site.position)),
                    'REF': "A",
                    'ALT': "G"
                },
                             index=[0]),
                ignore_index=True)

        snp_df = pd.DataFrame(d)
        return snp_df

    def mnms(self):
        """
        update .geno and .snp output to include simulated MNMs
        """

        snp_mnm = pd.DataFrame(
            columns=['ID', 'CHR', 'POS1', 'POS', 'REF', 'ALT'])
        #         geno_mat = self.geno()
        repeat_rows = []
        for index, snp in self.snp.iterrows():
            snp_mnm = snp_mnm.append(snp, ignore_index=True)

            random.seed(snp['POS'])
            if random.random() < self.mnm_frac:

                dist = random.randint(1, self.mnm_dist)

                mnm_snp = snp
                mnm_snp['POS'] = str(round(int(snp['POS']) + dist))
                mnm_snp['POS1'] = float(mnm_snp['POS']) / 10e6
                mnm_snp['ID'] = "1:" + mnm_snp['POS']
                snp_mnm = snp_mnm.append(mnm_snp, ignore_index=True)
                repeat_rows.append(2)
            else:
                repeat_rows.append(1)

        geno_mat_mnm = np.repeat(self.geno, repeats=repeat_rows, axis=0)

        mnm_prefix = self.prefix + "_mnm" + str(self.mnm_dist) + "-" + str(
            self.mnm_frac)

        return mnmEigData(
            geno=geno_mat_mnm,
            snp=snp_mnm,
            prefix=mnm_prefix,
            repeat_rows=repeat_rows)


class writeEigData:
    """
    functions for writing .snp .geno and .ind files
    """

    def __init__(self, eig_data, output_dir="./"):
        self.eig_data = eig_data
        self.prefix = eig_data.prefix
        self.output_dir = output_dir
        self.out_file = self.output_dir + self.prefix + ".snp"

    def dump(self):
        """
        run all 3 output functions: write_snp, write_geno, and write_ind
        """

        self.write_snp()
        self.write_geno()
        self.write_ind()

    def write_snp(self):
        """
        write .snp files
        """

        snp_out_file = self.output_dir + self.prefix + ".snp"
        util_log.debug("writing .snp file to %s", snp_out_file)
        self.eig_data.snp.to_csv(snp_out_file, index=False, header=False, sep="\t")

    def write_geno(self):
        """
        write .geno files
        """

        geno_afr = self.eig_data.geno[:, :100]
        geno_eur = self.eig_data.geno[:, 100:200]

        geno_pop = [geno_afr, geno_eur]

        for i, pop in enumerate(["afr", "eur"]):
            geno_out_file = self.output_dir + self.prefix + "_" + pop + ".geno"
            util_log.debug("writing .geno file to %s", geno_out_file)
            np.savetxt(geno_out_file, geno_pop[i], delimiter="", fmt='%i')

    def write_ind(self):
        """
        write .ind files
        """

        for pop in ["afr", "eur"]:
            # write out separate .ind files per population.
            # columns indicate sample ID, sex (set as 'U'), and label (set as 'ADMIXED')
            ind_file = self.output_dir + self.prefix + "_" + pop + ".ind"
            util_log.debug("writing .ind file to %s", ind_file)
            with open(ind_file, "w") as id_file:
                for sample_id in range(0, 100):
                    sample_name = self.prefix + "_" + pop + \
                        "_sample_" + str(sample_id)
                    print("\t".join([sample_name, "U", pop]), file=id_file)
