import os
import sys
import argparse
import warnings
import multiprocessing
import msprime
# import pandas as pd
from joblib import Parallel, delayed

import util

proj_dir = os.getcwd()
sys.path.insert(0, proj_dir + '/src/stdpopsim')

from stdpopsim import homo_sapiens, models


def main():

    parser = argparse.ArgumentParser()

    num_cores = multiprocessing.cpu_count()
    parser.add_argument(
        "-c",
        "--cpus",
        help="number of CPUs. Must be integer value between 1 \
                            and " + str(num_cores),
        nargs='?',
        type=int,
        choices=range(1, num_cores + 1),
        metavar='INT',
        default=num_cores)

    parser.add_argument(
        "-v", "--verbose", help="Enable verbose logging", action="store_true")

    parser.add_argument(
        "-m",
        "--mutation_rate",
        help="mutation rate in simulation",
        nargs='?',
        type=float,
        metavar='FLOAT',
        default=1.15e-8)

    parser.add_argument(
        "-r",
        "--recombination_rate",
        help="recombination rate in simulation",
        nargs='?',
        type=float,
        metavar='FLOAT',
        default=1e-8)

    parser.add_argument(
        "-f",
        "--mnm_frac",
        help="fraction of expected MNMs",
        nargs='?',
        type=float,
        metavar='FLOAT',
        default=0.015)

    parser.add_argument(
        "-d",
        "--mnm_dist",
        help="maximum distance between simulated MNMs",
        nargs='?',
        type=int,
        metavar='INT',
        default=100)

    parser.add_argument(
        "-l",
        "--length",
        help="length of each simulated haplotype",
        nargs='?',
        type=int,
        metavar='INT',
        default=50000)

    parser.add_argument(
        "-N",
        "--num_samples",
        help="number of samples per replicate",
        nargs='?',
        type=int,
        metavar='INT',
        default=100)

    parser.add_argument(
        "-R",
        "--replicates",
        help="number of replicates",
        nargs='?',
        type=int,
        metavar='INT',
        default=1000)

    parser.add_argument(
        "-i",
        "--replicate_ID",
        help="unique identifier when running simulation for specific replicate",
        nargs='?',
        type=int,
        metavar='INT',
        default=0)

    parser.add_argument(
        "-s",
        "--seed",
        help="set seed",
        nargs='?',
        type=int,
        metavar='INT',
        default=30)

    parser.add_argument(
        "-D",
        "--demographic_model",
        help="demographic model to simulate under",
        nargs='?',
        type=str,
        metavar='STR',
        default="GutenkunstThreePop")

    parser.add_argument(
        "-M",
        "--method",
        help="archaic ancestry inference method",
        nargs='?',
        type=str,
        metavar='method',
        default="archie")

    parser.add_argument(
        "-F",
        "--force",
        help="force rewrite of output files, even if they already exist",
        action="store_true")

    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    if args.verbose:
        loglev = 'DEBUG'
    else:
        loglev = 'INFO'
        warnings.filterwarnings("ignore", category=UserWarning)

    util.util_log.setLevel(loglev)
    log = util.get_logger("archanc", level=loglev)

    log.debug("Running with the following options:")
    for arg in vars(args):
        log.debug("%s : %s", arg, getattr(args, arg))

    # coalescent simulation parameters
    sample_size = args.num_samples
    length = args.length
    mu = args.mutation_rate
    rr = args.recombination_rate
    replicates = args.replicates

    if replicates == 1:
        seed = args.replicate_ID
    else:
        seed = args.seed

    proj_dir = os.getcwd()
    msprime_dir = proj_dir + "/output/msprime/"
    out_dir = msprime_dir + args.demographic_model + "/"
    archie_src_dir = proj_dir + "/src/ArchIE/"
    archie_out_dir = proj_dir + "/output/ArchIE/"

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if args.demographic_model == "GutenkunstThreePop":

        demo_model = homo_sapiens.GutenkunstThreePopOutOfAfrica()
        pop_samples = [msprime.Sample(0, 0)] * sample_size + \
                [msprime.Sample(1, 0)] * sample_size + \
                [msprime.Sample(2, 0)] * sample_size

        # model_dict = {"GutenkunstThreePop": demo_model_ts}

    elif args.demographic_model == "TennessenTwoPop":

        demo_model = homo_sapiens.TennessenTwoPopOutOfAfrica()
        pop_samples = [msprime.Sample(0, 0)] * sample_size + \
                [msprime.Sample(1, 0)] * sample_size

    elif args.demographic_model == "RagsdaleArchaic":

        demo_model = homo_sapiens.RagsdaleArchaic()
        pop_samples = [msprime.Sample(0, 0)] * sample_size + \
                [msprime.Sample(1, 0)] * sample_size

    demo_model_ts = msprime.simulate(
        # first 100 samples from AFR, next 100 from EUR
        samples=pop_samples,
        length=length,
        mutation_rate=mu,
        recombination_rate=rr,
        random_seed=seed,
        num_replicates=replicates,
        **demo_model.asdict())

    model_dict = {"model": demo_model_ts}

    #-------------------------------------------------------
    # define other models here and add to model_dict below
    #-------------------------------------------------------
    # e.g., modify demographic parameters to include archaic branches
    # GutenkunstThreePopArchaic_model = homo_sapiens.GutenkunstThreePopArchaic()

    # GutenkunstThreePopArchaic_ts = msprime.simulate(
    #     # first 100 samples from AFR, next 100 from EUR
    #     samples=[msprime.Sample(0, 0)]*sample_size + [msprime.Sample(1, 0)]*sample_size,
    #     length=length,
    #     mutation_rate=mu,
    #     recombination_rate=rr,
    #     random_seed=seed,
    #     num_replicates=replicates,
    #     **GutenkunstThreePopArchaic_model.asdict())

    # process simulated data and write files for use with different methods
    for model_label, model in model_dict.items():

        ts_list = {}
        for j, ts in enumerate(model):

            if args.replicate_ID == 0:
                prefix_nomnm = out_dir + args.demographic_model + "_rep" + str(
                    j + 1)
            else:
                prefix_nomnm = out_dir + args.demographic_model + "_rep" + str(
                    args.replicate_ID)
            prefix_mnm = prefix_nomnm + "_mnm" + str(
                args.mnm_dist) + "-" + str(args.mnm_frac)

            if args.method == "sprime":
                suffix = ".vcf"

            elif args.method == "archie":
                suffix = ".snp"

            if args.force or not os.path.isfile(
                    prefix_nomnm + suffix) or not os.path.isfile(prefix_mnm +
                                                                 suffix):
                log.info(
                    "Output files %s*%s are missing or the --force flag is enabled"
                    % (prefix_nomnm, suffix))
                ts_list[j] = list(ts.variants())
            else:
                log.debug(
                    "Output files %s*%s already exist and will not be overwritten"
                    % (prefix_nomnm, suffix))

        # parallelize if running multiple replicates
        if args.cpus > 1 and replicates > 1:
            Parallel(n_jobs=args.cpus) \
                (delayed(util.process_ts)(ts, args.demographic_model, j+1, args.mnm_frac, args.mnm_dist, args.method, out_dir) \
                for j, ts in ts_list.items())

        # run as single instance if returning only a single replicate
        # (for use with the cluster)
        elif replicates == 1:
            for j, ts in ts_list.items():
                util.process_ts(ts, args.demographic_model, args.replicate_ID,
                                args.mnm_frac, args.mnm_dist, args.method,
                                out_dir)

        # for j, ts in enumerate(model):
        #     util.process_ts(ts.variants(), model_label, j, args.mnm_frac, args.mnm_dist, args.method, out_dir)


#             if args.method == "archie":
#                 for data in [eig_data, eig_data_mnms]:
#                     for pop in ["afr", "eur"]:

#                         if pop == "afr":
#                             ref_pop = "eur"
#                         else:
#                             ref_pop = "afr"

#                         prefix = data.prefix
#                         stats_pop_cmd = "python " + archie_src_dir + "data/calc_stats_window_data.py" + \
#                             " -s " + msprime_dir + prefix + ".snp" + \
#                             " -i " + msprime_dir + prefix + "_" + pop + ".ind" + \
#                             " -a " + msprime_dir + prefix + "_" + pop + ".geno" + \
#                             " -r " + msprime_dir + prefix + "_" + ref_pop + ".geno" + \
#                             " -c 1 -b 0 -e 50000 -w 50000 -z 50000 " + \
#                             " > " + archie_out_dir + prefix + "_" + pop + ".txt"
#                         print(stats_pop_cmd + "\n")
# #                         os.system(stats_pop_cmd)

main()
