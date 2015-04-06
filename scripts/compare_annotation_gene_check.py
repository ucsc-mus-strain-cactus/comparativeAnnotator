#!/usr/bin/env python
import sys, os
fields = (("UtrGap", "utrGap"), ("CdsGap", "cdsGap"), ("BeginStart", "noStart"), ("EndStop", "noStop"), 
          ("CdsUnknownSplice", "unknownCdsSplice"), )

def main():
    mark_records = [x for x in open(sys.argv[1])]
    my_records = [x for x in open(sys.argv[2])]
    with open(sys.argv[3], "w") as outf:
        outf.write("Field\tMyValue\tMarkValue\n")
        for my_field, mark_field in fields:
            # stupid
            if my_field == "CdsGap":
                my_this = len([x for x in my_records if my_field in x or "CdsMult3Gap" in x])
            else:
                my_this = len([x for x in my_records if my_field in x])
            mark_this = len([x for x in mark_records if mark_field in x])
            outf.write("{}/{}\t{}\t{}\n".format(my_field, mark_field, my_this, mark_this))

    # hacked
    d = os.path.dirname(sys.argv[3])
    for my_field, mark_field in fields:
        with open(os.path.join(d, "mark_diff_my_" + my_field + ".bed"), "w") as outf:
            outf.write('track name="mark_diff_my_{}"\n'.format(my_field))
            my_this = {tuple(x.split()[:3] + [x.split()[3].split("/")[-1]]) for x in my_records if my_field in x}
            mark_this = {tuple(x.split()[:3] + [x.split()[3].split("/")[-1]]) for x in mark_records if mark_field in x}
            q = list(mark_this.difference(my_this))
            for x in q:
                outf.write("\t".join(x)+"\n")
        with open(os.path.join(d, "my_diff_mark_" + my_field + ".bed"), "w") as outf:
            outf.write('track name="my_diff_mark_{}"\n'.format(my_field))
            my_this = {tuple(x.split()[:3] + [x.split()[3].split("/")[-1]]) for x in my_records if my_field in x}
            mark_this = {tuple(x.split()[:3] + [x.split()[3].split("/")[-1]]) for x in mark_records if mark_field in x}
            q = list(my_this.difference(mark_this))
            for x in q:
                outf.write("\t".join(x)+"\n")        


if __name__ == '__main__':
    main()