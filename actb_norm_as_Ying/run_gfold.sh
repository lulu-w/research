gfold diff -norm NO -s1 sen_untreated -s2 rej_untreated -suf .read_cnt -o sen_untreatedVSrej_untreated.diff; 

gfold diff -norm NO -s1 sr_untreated -s2 sen_untreated -suf .read_cnt -o sr_untreatedVSsen_untreated.diff; gfold diff -norm NO -s1 sr_untreated -s2 rej_untreated -suf .read_cnt -o sr_untreatedVSrej_untreated.diff;