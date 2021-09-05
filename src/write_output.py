import multiprocessing as mp
from tqdm import tqdm
from transform_data import *

def write_ouput(de):
    if de.entropy:
        de.MOTUoutfile = str(de.MOTUoutfile + '_Adcorr')

    print('writing output_info')

    de.output_info = pd.DataFrame.from_dict(de.output_info)
    de.output_info.to_csv(str(de.MOTUoutfile + '_denoising_info.csv'), index=False)

    if (de.output_type == 'ratio') or (de.output_type == 'all'):
        mothers_ratio = de.output_info.mother_ratio.unique()[1:]
    if (de.output_type == 'd') or (de.output_type == 'all'):
        mothers_d = de.output_info.mother_d.unique()[1:]
    if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
        mothers_ratio_d = de.output_info.mother_xavier_criteria.unique()[1:]

    del de.output_info

    de.data_initial = de.data_initial.set_index(de.data_initial.loc[:, 'id'])

    if (de.output_type == 'ratio') or (de.output_type == 'all'):
        de.good_mothers = de.data_initial.loc[de.good_seq][de.first_col_names + de.abund_col_names + [de.seq]]
        print('writing output_ratio')
        # writing ratio
        if de.cores > 1:
            pool = mp.Pool(de.cores)
            [row] = zip(*pool.map(de.write_output_ratio, [mother for mother in mothers_ratio]))
            pool.close()
            del pool
            denoised_ratio = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            de.good_mothers = de.good_mothers.drop(index=mothers_ratio)
            denoised_ratio = denoised_ratio.append(de.good_mothers, ignore_index=True)
            denoised_ratio = denoised_ratio.sort_values([de.count], axis=0, ascending=False)
        else:
            denoised_ratio = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            for mother in tqdm(mothers_ratio):
                row = [
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[
                        0] +
                    list(de.data_initial.loc[
                             list(pd.Series(de.denoised_ratio_output) == mother), de.abund_col_names].sum(0)) +
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
                row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
                denoised_ratio = denoised_ratio.append(row, ignore_index=True)
                de.good_mothers = de.good_mothers.drop(index=mother)
            denoised_ratio = denoised_ratio.append(de.good_mothers, ignore_index=True)
            denoised_ratio = denoised_ratio.sort_values([de.count], axis=0, ascending=False)
        if 'row' in locals():
            del row
        del de.good_mothers, mothers_ratio, de.denoised_ratio_output

        if de.fasta_output:
            denoised_ratio = denoised_ratio.to_dict(orient='index')
            ofile = open(str(de.MOTUoutfile + '_denoised_ratio.fasta'), "w")
            for i in tqdm(range(len(denoised_ratio))):
                ofile.write(">" + denoised_ratio[i]['id'] + ';size=' + str(denoised_ratio[i][de.count]) +
                            ";\n" + denoised_ratio[i][de.seq].upper() + "\n")
            # do not forget to close it
            ofile.close()
        else:
            denoised_ratio.to_csv(str(de.MOTUoutfile + '_denoised_ratio.csv'), index=False)

        del denoised_ratio

    if (de.output_type == 'd') or (de.output_type == 'all'):
        de.good_mothers = de.data_initial.loc[de.good_seq][de.first_col_names + de.abund_col_names + [de.seq]]
        print('writing output_d')
        # writing d
        if de.cores > 1:
            pool = mp.Pool(de.cores)
            [row] = zip(*pool.map(de.write_output_d, [mother for mother in mothers_d]))
            pool.close()
            del pool
            denoised_d = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            de.good_mothers = de.good_mothers.drop(index=mothers_d)
            denoised_d = denoised_d.append(de.good_mothers, ignore_index=True)
            denoised_d = denoised_d.sort_values([de.count], axis=0, ascending=False)
        else:
            denoised_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            for mother in tqdm(mothers_d):
                row = [
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[
                        0] +
                    list(de.data_initial.loc[list(pd.Series(de.denoised_d_output) == mother), de.abund_col_names].sum(
                        0)) +
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
                row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
                denoised_d = denoised_d.append(row, ignore_index=True)
                de.good_mothers = de.good_mothers.drop(index=mother)
            denoised_d = denoised_d.append(de.good_mothers, ignore_index=True)
            denoised_d = denoised_d.sort_values([de.count], axis=0, ascending=False)
        if 'row' in locals():
            del row
        del de.good_mothers, mothers_d, de.denoised_d_output

        if de.fasta_output:
            denoised_d = denoised_d.to_dict(orient='index')
            ofile = open(str(de.MOTUoutfile + '_denoised_d.fasta'), "w")
            for i in tqdm(range(len(denoised_d))):
                ofile.write(">" + denoised_d[i]['id'] + ';size=' + str(denoised_d[i][de.count]) + ";\n" +
                            denoised_d[i][de.seq].upper() + "\n")
            # do not forget to close it
            ofile.close()
        else:
            denoised_d.to_csv(str(de.MOTUoutfile + '_denoised_d.csv'), index=False)

        del denoised_d

    if (de.output_type == 'ratio_d') or (de.output_type == 'all'):
        de.good_mothers = de.data_initial.loc[de.good_seq][de.first_col_names + de.abund_col_names + [de.seq]]
        print('writing output_ratio_d')
        # writing ratio_d
        if de.cores > 1:
            pool = mp.Pool(de.cores)
            [row] = zip(*pool.map(de.write_output_ratio_d, [mother for mother in mothers_ratio_d]))
            pool.close()
            del pool
            denoised_ratio_d = pd.DataFrame(row, columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            de.good_mothers = de.good_mothers.drop(index=mothers_ratio_d)
            denoised_ratio_d = denoised_ratio_d.append(de.good_mothers, ignore_index=True)
            denoised_ratio_d = denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
        else:
            denoised_ratio_d = pd.DataFrame(columns=[de.first_col_names + de.abund_col_names + [de.seq]][0])
            for mother in tqdm(mothers_ratio_d):
                row = [
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.first_col_names].values.tolist()[
                        0] +
                    list(de.data_initial.loc[
                        list(pd.Series(de.denoised_ratio_d_output) == mother), de.abund_col_names].sum(
                        0)) +
                    de.good_mothers[list(de.good_mothers.loc[:, 'id'] == mother)][de.seq].values.tolist()]
                row = pd.Series(row[0], index=[de.first_col_names + de.abund_col_names + [de.seq]][0])
                denoised_ratio_d = denoised_ratio_d.append(row, ignore_index=True)
                de.good_mothers = de.good_mothers.drop(index=mother)
            denoised_ratio_d = denoised_ratio_d.append(de.good_mothers, ignore_index=True)
            denoised_ratio_d = denoised_ratio_d.sort_values([de.count], axis=0, ascending=False)
        if 'row' in locals():
            del row
        del de.good_mothers, mothers_ratio_d, de.denoised_ratio_d_output

        if de.fasta_output:
            denoised_ratio_d = denoised_ratio_d.to_dict(orient='index')
            ofile = open(str(de.MOTUoutfile + '_denoised_ratio_d.fasta'), "w")
            for i in tqdm(range(len(denoised_ratio_d))):
                ofile.write(">" + denoised_ratio_d[i]['id'] + ';size=' + str(denoised_ratio_d[i][de.count]) + ";\n" +
                            denoised_ratio_d[i][de.seq].upper() + "\n")
            # do not forget to close it
            ofile.close()
        else:
            denoised_ratio_d.to_csv(str(de.MOTUoutfile + '_denoised_ratio_d.csv'), index=False)

        del denoised_ratio_d

