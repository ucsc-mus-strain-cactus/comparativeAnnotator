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

if __name__ == '__main__':
    main()