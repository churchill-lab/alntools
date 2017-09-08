/**
 * An ordered array of chromosome sizes. This data structure was derived from the data in
 * ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz
 * @type Array.<GenoInterval>
 */
var hg19ChrSizes = [
    {
        "chr": "1",
        "size": 249250621,
        "startPos": 0
    },
    {
        "chr": "2",
        "size": 243199373,
        "startPos": 0
    },
    {
        "chr": "3",
        "size": 198022430,
        "startPos": 0
    },
    {
        "chr": "4",
        "size": 191154276,
        "startPos": 0
    },
    {
        "chr": "5",
        "size": 180915260,
        "startPos": 0
    },
    {
        "chr": "6",
        "size": 171115067,
        "startPos": 0
    },
    {
        "chr": "7",
        "size": 159138663,
        "startPos": 0
    },
    {
        "chr": "8",
        "size": 146364022,
        "startPos": 0
    },
    {
        "chr": "9",
        "size": 141213431,
        "startPos": 0
    },
    {
        "chr": "10",
        "size": 135534747,
        "startPos": 0
    },
    {
        "chr": "11",
        "size": 135006516,
        "startPos": 0
    },
    {
        "chr": "12",
        "size": 133851895,
        "startPos": 0
    },
    {
        "chr": "13",
        "size": 115169878,
        "startPos": 0
    },
    {
        "chr": "14",
        "size": 107349540,
        "startPos": 0
    },
    {
        "chr": "15",
        "size": 102531392,
        "startPos": 0
    },
    {
        "chr": "16",
        "size": 90354753,
        "startPos": 0
    },
    {
        "chr": "17",
        "size": 81195210,
        "startPos": 0
    },
    {
        "chr": "18",
        "size": 78077248,
        "startPos": 0
    },
    {
        "chr": "19",
        "size": 59128983,
        "startPos": 0
    },
    {
        "chr": "20",
        "size": 63025520,
        "startPos": 0
    },
    {
        "chr": "21",
        "size": 48129895,
        "startPos": 0
    },
    {
        "chr": "22",
        "size": 51304566,
        "startPos": 0
    },
    {
        "chr": "X",
        "size": 155270560,
        "startPos": 0
    },
    {
        "chr": "Y",
        "size": 59373566,
        "startPos": 0
    },
    {
        "chr": "M",
        "size": 16571,
        "startPos": 0
    }
];

/**
 * An array of cytoband intervals. Note that in addition to the GenoInterval properties each array element has a
 * cytoBandType string property (eg. "gpos66"). This data structure was derived from the data in
 * ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
 * @type Array.<GenoInterval>
 */
var hg19CytoBands = [
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 0
    },
    {
        "chr": "1",
        "cytoBandType": "gpos25",
        "size": 3100000,
        "startPos": 2300000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 5400000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos25",
        "size": 2000000,
        "startPos": 7200000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 3500000,
        "startPos": 9200000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 3500000,
        "startPos": 12700000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 4200000,
        "startPos": 16200000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos25",
        "size": 3500000,
        "startPos": 20400000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 4100000,
        "startPos": 23900000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos25",
        "size": 2200000,
        "startPos": 28000000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 30200000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos25",
        "size": 2200000,
        "startPos": 32400000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 5500000,
        "startPos": 34600000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos25",
        "size": 4000000,
        "startPos": 40100000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 44100000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos75",
        "size": 3900000,
        "startPos": 46800000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 5400000,
        "startPos": 50700000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 2900000,
        "startPos": 56100000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 59000000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 7600000,
        "startPos": 61300000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 800000,
        "startPos": 68900000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos100",
        "size": 15200000,
        "startPos": 69700000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 3500000,
        "startPos": 84900000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos75",
        "size": 3600000,
        "startPos": 88400000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 92000000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos75",
        "size": 5000000,
        "startPos": 94700000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 2500000,
        "startPos": 99700000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos100",
        "size": 5000000,
        "startPos": 102200000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 4600000,
        "startPos": 107200000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 4300000,
        "startPos": 111800000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 116100000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 2800000,
        "startPos": 117800000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 900000,
        "startPos": 120600000
    },
    {
        "chr": "1",
        "cytoBandType": "acen",
        "size": 3500000,
        "startPos": 121500000
    },
    {
        "chr": "1",
        "cytoBandType": "acen",
        "size": 3900000,
        "startPos": 125000000
    },
    {
        "chr": "1",
        "cytoBandType": "gvar",
        "size": 13700000,
        "startPos": 128900000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 142600000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 3300000,
        "startPos": 147000000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 4700000,
        "startPos": 150300000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 1500000,
        "startPos": 155000000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 156500000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 1400000,
        "startPos": 159100000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 5000000,
        "startPos": 160500000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 1700000,
        "startPos": 165500000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 3700000,
        "startPos": 167200000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos75",
        "size": 2000000,
        "startPos": 170900000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 3100000,
        "startPos": 172900000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 4300000,
        "startPos": 176000000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 5500000,
        "startPos": 180300000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos100",
        "size": 5000000,
        "startPos": 185800000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 3000000,
        "startPos": 190800000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos100",
        "size": 4900000,
        "startPos": 193800000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 8500000,
        "startPos": 198700000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos25",
        "size": 4300000,
        "startPos": 207200000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 3000000,
        "startPos": 211500000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos100",
        "size": 9600000,
        "startPos": 214500000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 500000,
        "startPos": 224100000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos25",
        "size": 2400000,
        "startPos": 224600000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 3700000,
        "startPos": 227000000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos50",
        "size": 4000000,
        "startPos": 230700000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 234700000
    },
    {
        "chr": "1",
        "cytoBandType": "gpos75",
        "size": 7100000,
        "startPos": 236600000
    },
    {
        "chr": "1",
        "cytoBandType": "gneg",
        "size": 5550621,
        "startPos": 243700000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 0
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 2700000,
        "startPos": 4400000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 5100000,
        "startPos": 7100000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos75",
        "size": 4500000,
        "startPos": 12200000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 2500000,
        "startPos": 16700000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos75",
        "size": 4800000,
        "startPos": 19200000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 3900000,
        "startPos": 24000000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos25",
        "size": 2100000,
        "startPos": 27900000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 30000000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos75",
        "size": 4500000,
        "startPos": 32100000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 36600000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 3200000,
        "startPos": 38600000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 6000000,
        "startPos": 41800000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos100",
        "size": 5100000,
        "startPos": 47800000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 52900000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos100",
        "size": 6300000,
        "startPos": 55000000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 61300000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 4500000,
        "startPos": 64100000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 68600000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 2000000,
        "startPos": 71500000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 1500000,
        "startPos": 73500000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos100",
        "size": 8300000,
        "startPos": 75000000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 7200000,
        "startPos": 83300000
    },
    {
        "chr": "2",
        "cytoBandType": "acen",
        "size": 2800000,
        "startPos": 90500000
    },
    {
        "chr": "2",
        "cytoBandType": "acen",
        "size": 3500000,
        "startPos": 93300000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 5900000,
        "startPos": 96800000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 3300000,
        "startPos": 102700000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 1500000,
        "startPos": 106000000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos25",
        "size": 2700000,
        "startPos": 107500000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 4200000,
        "startPos": 110200000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 4400000,
        "startPos": 114400000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 3600000,
        "startPos": 118800000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 7500000,
        "startPos": 122400000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 129900000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos25",
        "size": 2600000,
        "startPos": 132500000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 135100000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos100",
        "size": 5400000,
        "startPos": 136800000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 142200000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos100",
        "size": 4600000,
        "startPos": 144100000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 1200000,
        "startPos": 148700000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos25",
        "size": 600000,
        "startPos": 149900000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 150500000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos75",
        "size": 4900000,
        "startPos": 154900000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 3900000,
        "startPos": 159800000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos75",
        "size": 6000000,
        "startPos": 163700000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 8300000,
        "startPos": 169700000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 2600000,
        "startPos": 178000000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 2400000,
        "startPos": 180600000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos75",
        "size": 6400000,
        "startPos": 183000000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 2500000,
        "startPos": 189400000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos75",
        "size": 5500000,
        "startPos": 191900000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 5900000,
        "startPos": 197400000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 1600000,
        "startPos": 203300000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 4100000,
        "startPos": 204900000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos100",
        "size": 6300000,
        "startPos": 209000000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 6200000,
        "startPos": 215300000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos75",
        "size": 3700000,
        "startPos": 221500000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 900000,
        "startPos": 225200000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos100",
        "size": 4900000,
        "startPos": 226100000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 4600000,
        "startPos": 231000000
    },
    {
        "chr": "2",
        "cytoBandType": "gpos50",
        "size": 1700000,
        "startPos": 235600000
    },
    {
        "chr": "2",
        "cytoBandType": "gneg",
        "size": 5899373,
        "startPos": 237300000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos50",
        "size": 2800000,
        "startPos": 0
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1200000,
        "startPos": 2800000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos50",
        "size": 4700000,
        "startPos": 4000000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 3100000,
        "startPos": 8700000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos25",
        "size": 1500000,
        "startPos": 11800000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 3100000,
        "startPos": 13300000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos100",
        "size": 7500000,
        "startPos": 16400000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 2500000,
        "startPos": 23900000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 4500000,
        "startPos": 26400000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1200000,
        "startPos": 30900000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos50",
        "size": 4400000,
        "startPos": 32100000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 36500000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 4300000,
        "startPos": 39400000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 400000,
        "startPos": 43700000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos50",
        "size": 100000,
        "startPos": 44100000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 6400000,
        "startPos": 44200000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos25",
        "size": 1700000,
        "startPos": 50600000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 52300000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos50",
        "size": 4200000,
        "startPos": 54400000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 5100000,
        "startPos": 58600000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos50",
        "size": 6100000,
        "startPos": 63700000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 69800000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 5600000,
        "startPos": 74200000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 3700000,
        "startPos": 79800000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 3700000,
        "startPos": 83500000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 700000,
        "startPos": 87200000
    },
    {
        "chr": "3",
        "cytoBandType": "acen",
        "size": 3100000,
        "startPos": 87900000
    },
    {
        "chr": "3",
        "cytoBandType": "acen",
        "size": 2900000,
        "startPos": 91000000
    },
    {
        "chr": "3",
        "cytoBandType": "gvar",
        "size": 4400000,
        "startPos": 93900000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 98300000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos25",
        "size": 900000,
        "startPos": 100000000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 100900000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 3400000,
        "startPos": 102800000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 106200000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos50",
        "size": 3400000,
        "startPos": 107900000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 111300000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 3800000,
        "startPos": 113500000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 117300000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 2900000,
        "startPos": 119000000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 121900000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos25",
        "size": 2000000,
        "startPos": 123800000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 125800000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos25",
        "size": 4500000,
        "startPos": 129200000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 133700000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos25",
        "size": 3000000,
        "startPos": 135700000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 4100000,
        "startPos": 138700000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos100",
        "size": 6100000,
        "startPos": 142800000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 3200000,
        "startPos": 148900000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos50",
        "size": 2900000,
        "startPos": 152100000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 155000000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos50",
        "size": 2000000,
        "startPos": 157000000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 159000000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos100",
        "size": 6900000,
        "startPos": 160700000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 167600000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 4800000,
        "startPos": 170900000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 175700000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 3700000,
        "startPos": 179000000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 182700000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos25",
        "size": 1500000,
        "startPos": 184500000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 186000000
    },
    {
        "chr": "3",
        "cytoBandType": "gpos75",
        "size": 4400000,
        "startPos": 187900000
    },
    {
        "chr": "3",
        "cytoBandType": "gneg",
        "size": 5722430,
        "startPos": 192300000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 4500000,
        "startPos": 0
    },
    {
        "chr": "4",
        "cytoBandType": "gpos25",
        "size": 1500000,
        "startPos": 4500000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 5300000,
        "startPos": 6000000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos50",
        "size": 3900000,
        "startPos": 11300000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 15200000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos75",
        "size": 3500000,
        "startPos": 17800000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 6400000,
        "startPos": 21300000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos100",
        "size": 8100000,
        "startPos": 27700000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 5400000,
        "startPos": 35800000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos50",
        "size": 3400000,
        "startPos": 41200000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 3600000,
        "startPos": 44600000
    },
    {
        "chr": "4",
        "cytoBandType": "acen",
        "size": 2200000,
        "startPos": 48200000
    },
    {
        "chr": "4",
        "cytoBandType": "acen",
        "size": 2300000,
        "startPos": 50400000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 6800000,
        "startPos": 52700000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos100",
        "size": 7100000,
        "startPos": 59500000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 3900000,
        "startPos": 66600000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos75",
        "size": 5800000,
        "startPos": 70500000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 76300000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos50",
        "size": 3500000,
        "startPos": 78900000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 82400000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos25",
        "size": 2800000,
        "startPos": 84100000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 1100000,
        "startPos": 86900000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos75",
        "size": 5700000,
        "startPos": 88000000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 1400000,
        "startPos": 93700000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos75",
        "size": 3700000,
        "startPos": 95100000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 98800000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos50",
        "size": 6600000,
        "startPos": 101100000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 6400000,
        "startPos": 107700000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos75",
        "size": 6700000,
        "startPos": 114100000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 3000000,
        "startPos": 120800000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos50",
        "size": 5000000,
        "startPos": 123800000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 128800000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos100",
        "size": 8400000,
        "startPos": 131100000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 139500000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos25",
        "size": 5300000,
        "startPos": 141500000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 146800000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos25",
        "size": 2600000,
        "startPos": 148500000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 4500000,
        "startPos": 151100000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos100",
        "size": 6200000,
        "startPos": 155600000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 161800000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos100",
        "size": 5600000,
        "startPos": 164500000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 170100000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos75",
        "size": 4400000,
        "startPos": 171900000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 1200000,
        "startPos": 176300000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos100",
        "size": 5700000,
        "startPos": 177500000
    },
    {
        "chr": "4",
        "cytoBandType": "gneg",
        "size": 3900000,
        "startPos": 183200000
    },
    {
        "chr": "4",
        "cytoBandType": "gpos25",
        "size": 4054276,
        "startPos": 187100000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 4500000,
        "startPos": 0
    },
    {
        "chr": "5",
        "cytoBandType": "gpos25",
        "size": 1800000,
        "startPos": 4500000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 3500000,
        "startPos": 6300000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos50",
        "size": 5200000,
        "startPos": 9800000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 15000000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos100",
        "size": 4900000,
        "startPos": 18400000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 1300000,
        "startPos": 23300000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos100",
        "size": 4300000,
        "startPos": 24600000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 4900000,
        "startPos": 28900000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos25",
        "size": 4600000,
        "startPos": 33800000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 4100000,
        "startPos": 38400000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos50",
        "size": 3600000,
        "startPos": 42500000
    },
    {
        "chr": "5",
        "cytoBandType": "acen",
        "size": 2300000,
        "startPos": 46100000
    },
    {
        "chr": "5",
        "cytoBandType": "acen",
        "size": 2300000,
        "startPos": 48400000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 8200000,
        "startPos": 50700000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos75",
        "size": 4000000,
        "startPos": 58900000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 300000,
        "startPos": 62900000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos75",
        "size": 3500000,
        "startPos": 63200000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 66700000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos50",
        "size": 4900000,
        "startPos": 68400000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 3600000,
        "startPos": 73300000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos50",
        "size": 4500000,
        "startPos": 76900000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 1400000,
        "startPos": 81400000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos100",
        "size": 9500000,
        "startPos": 82800000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 5900000,
        "startPos": 92300000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos100",
        "size": 4600000,
        "startPos": 98200000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 102800000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos100",
        "size": 5100000,
        "startPos": 104500000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 109600000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos50",
        "size": 1600000,
        "startPos": 111500000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 113100000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos100",
        "size": 6200000,
        "startPos": 115200000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 5900000,
        "startPos": 121400000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos100",
        "size": 3300000,
        "startPos": 127300000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 5600000,
        "startPos": 130600000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos25",
        "size": 3300000,
        "startPos": 136200000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 5000000,
        "startPos": 139500000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos75",
        "size": 5300000,
        "startPos": 144500000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 149800000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos50",
        "size": 3000000,
        "startPos": 152700000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 4200000,
        "startPos": 155700000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos100",
        "size": 8600000,
        "startPos": 159900000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 4300000,
        "startPos": 168500000
    },
    {
        "chr": "5",
        "cytoBandType": "gpos25",
        "size": 3800000,
        "startPos": 172800000
    },
    {
        "chr": "5",
        "cytoBandType": "gneg",
        "size": 4315260,
        "startPos": 176600000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 0
    },
    {
        "chr": "6",
        "cytoBandType": "gpos25",
        "size": 1900000,
        "startPos": 2300000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 4200000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos50",
        "size": 3500000,
        "startPos": 7100000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 1000000,
        "startPos": 10600000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos25",
        "size": 1800000,
        "startPos": 11600000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 13400000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos75",
        "size": 10000000,
        "startPos": 15200000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 25200000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos50",
        "size": 3400000,
        "startPos": 27000000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 30400000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos25",
        "size": 1400000,
        "startPos": 32100000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 3100000,
        "startPos": 33500000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos25",
        "size": 3900000,
        "startPos": 36600000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 5700000,
        "startPos": 40500000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos100",
        "size": 5600000,
        "startPos": 46200000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 1100000,
        "startPos": 51800000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos100",
        "size": 4100000,
        "startPos": 52900000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 57000000
    },
    {
        "chr": "6",
        "cytoBandType": "acen",
        "size": 2300000,
        "startPos": 58700000
    },
    {
        "chr": "6",
        "cytoBandType": "acen",
        "size": 2300000,
        "startPos": 61000000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 100000,
        "startPos": 63300000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos100",
        "size": 6600000,
        "startPos": 63400000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 5900000,
        "startPos": 70000000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos50",
        "size": 8000000,
        "startPos": 75900000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 1000000,
        "startPos": 83900000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos50",
        "size": 3100000,
        "startPos": 84900000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 5100000,
        "startPos": 88000000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos100",
        "size": 6400000,
        "startPos": 93100000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 1100000,
        "startPos": 99500000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos100",
        "size": 4900000,
        "startPos": 100600000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 9100000,
        "startPos": 105500000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos75",
        "size": 3700000,
        "startPos": 114600000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 200000,
        "startPos": 118300000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos100",
        "size": 7600000,
        "startPos": 118500000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 1000000,
        "startPos": 126100000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos75",
        "size": 3200000,
        "startPos": 127100000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 900000,
        "startPos": 130300000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos50",
        "size": 4000000,
        "startPos": 131200000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 3800000,
        "startPos": 135200000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos75",
        "size": 3800000,
        "startPos": 139000000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 142800000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos75",
        "size": 3400000,
        "startPos": 145600000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 3500000,
        "startPos": 149000000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos50",
        "size": 3000000,
        "startPos": 152500000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 5500000,
        "startPos": 155500000
    },
    {
        "chr": "6",
        "cytoBandType": "gpos50",
        "size": 3500000,
        "startPos": 161000000
    },
    {
        "chr": "6",
        "cytoBandType": "gneg",
        "size": 6615067,
        "startPos": 164500000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 0
    },
    {
        "chr": "7",
        "cytoBandType": "gpos25",
        "size": 1700000,
        "startPos": 2800000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 4500000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos100",
        "size": 6500000,
        "startPos": 7300000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 13800000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos100",
        "size": 4400000,
        "startPos": 16500000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 4600000,
        "startPos": 20900000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos50",
        "size": 2500000,
        "startPos": 25500000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 800000,
        "startPos": 28000000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 6200000,
        "startPos": 28800000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 35000000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 6100000,
        "startPos": 37200000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 43300000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 3600000,
        "startPos": 45400000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 1500000,
        "startPos": 49000000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 3500000,
        "startPos": 50500000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 4000000,
        "startPos": 54000000
    },
    {
        "chr": "7",
        "cytoBandType": "acen",
        "size": 1900000,
        "startPos": 58000000
    },
    {
        "chr": "7",
        "cytoBandType": "acen",
        "size": 1800000,
        "startPos": 59900000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 5300000,
        "startPos": 61700000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos50",
        "size": 5200000,
        "startPos": 67000000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 5300000,
        "startPos": 72200000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos100",
        "size": 8900000,
        "startPos": 77500000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 86400000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 2900000,
        "startPos": 88200000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 91100000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 5200000,
        "startPos": 92800000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 5800000,
        "startPos": 98000000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos50",
        "size": 700000,
        "startPos": 103800000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 104500000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 7200000,
        "startPos": 107400000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 114600000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 3700000,
        "startPos": 117400000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 121100000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 3300000,
        "startPos": 123800000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 127100000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos25",
        "size": 1200000,
        "startPos": 129200000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 130400000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos50",
        "size": 5600000,
        "startPos": 132600000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 4900000,
        "startPos": 138200000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos75",
        "size": 4800000,
        "startPos": 143100000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 4700000,
        "startPos": 147900000
    },
    {
        "chr": "7",
        "cytoBandType": "gpos25",
        "size": 2500000,
        "startPos": 152600000
    },
    {
        "chr": "7",
        "cytoBandType": "gneg",
        "size": 4038663,
        "startPos": 155100000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 0
    },
    {
        "chr": "8",
        "cytoBandType": "gpos75",
        "size": 4000000,
        "startPos": 2200000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 6500000,
        "startPos": 6200000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos100",
        "size": 6300000,
        "startPos": 12700000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 4300000,
        "startPos": 19000000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos50",
        "size": 4100000,
        "startPos": 23300000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 1400000,
        "startPos": 27400000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos75",
        "size": 7700000,
        "startPos": 28800000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 36500000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos25",
        "size": 1400000,
        "startPos": 38300000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 39700000
    },
    {
        "chr": "8",
        "cytoBandType": "acen",
        "size": 2500000,
        "startPos": 43100000
    },
    {
        "chr": "8",
        "cytoBandType": "acen",
        "size": 2500000,
        "startPos": 45600000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 4100000,
        "startPos": 48100000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos75",
        "size": 400000,
        "startPos": 52200000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 52600000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos50",
        "size": 6100000,
        "startPos": 55500000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 600000,
        "startPos": 61600000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos50",
        "size": 3800000,
        "startPos": 62200000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 66000000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos50",
        "size": 2500000,
        "startPos": 68000000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 70500000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos100",
        "size": 4400000,
        "startPos": 73900000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 78300000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos75",
        "size": 4500000,
        "startPos": 80100000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 84600000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos100",
        "size": 6400000,
        "startPos": 86900000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 5700000,
        "startPos": 93300000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos25",
        "size": 2600000,
        "startPos": 99000000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 4600000,
        "startPos": 101600000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos75",
        "size": 4300000,
        "startPos": 106200000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 1600000,
        "startPos": 110500000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos100",
        "size": 5600000,
        "startPos": 112100000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 1500000,
        "startPos": 117700000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos50",
        "size": 3300000,
        "startPos": 119200000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 4800000,
        "startPos": 122500000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos50",
        "size": 4200000,
        "startPos": 127300000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 4900000,
        "startPos": 131500000
    },
    {
        "chr": "8",
        "cytoBandType": "gpos75",
        "size": 3500000,
        "startPos": 136400000
    },
    {
        "chr": "8",
        "cytoBandType": "gneg",
        "size": 6464022,
        "startPos": 139900000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 0
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 2400000,
        "startPos": 2200000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 4600000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos75",
        "size": 5200000,
        "startPos": 9000000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 2400000,
        "startPos": 14200000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 1900000,
        "startPos": 16600000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 1400000,
        "startPos": 18500000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos100",
        "size": 5700000,
        "startPos": 19900000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 2400000,
        "startPos": 25600000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos100",
        "size": 5200000,
        "startPos": 28000000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 3100000,
        "startPos": 33200000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 2100000,
        "startPos": 36300000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 38400000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos50",
        "size": 2600000,
        "startPos": 41000000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 3700000,
        "startPos": 43600000
    },
    {
        "chr": "9",
        "cytoBandType": "acen",
        "size": 1700000,
        "startPos": 47300000
    },
    {
        "chr": "9",
        "cytoBandType": "acen",
        "size": 1700000,
        "startPos": 49000000
    },
    {
        "chr": "9",
        "cytoBandType": "gvar",
        "size": 15200000,
        "startPos": 50700000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 65900000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 3500000,
        "startPos": 68700000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 72200000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos50",
        "size": 5200000,
        "startPos": 74000000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 79200000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos50",
        "size": 3000000,
        "startPos": 81100000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 84100000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos50",
        "size": 3500000,
        "startPos": 86900000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 1400000,
        "startPos": 90400000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 2100000,
        "startPos": 91800000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 93900000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 2700000,
        "startPos": 96600000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 99300000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos100",
        "size": 5600000,
        "startPos": 102600000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 3100000,
        "startPos": 108200000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 3600000,
        "startPos": 111300000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 114900000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos75",
        "size": 4800000,
        "startPos": 117700000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 122500000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 4500000,
        "startPos": 125800000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 3200000,
        "startPos": 130300000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 500000,
        "startPos": 133500000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 134000000
    },
    {
        "chr": "9",
        "cytoBandType": "gpos25",
        "size": 1500000,
        "startPos": 135900000
    },
    {
        "chr": "9",
        "cytoBandType": "gneg",
        "size": 3813431,
        "startPos": 137400000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 3000000,
        "startPos": 0
    },
    {
        "chr": "10",
        "cytoBandType": "gpos25",
        "size": 800000,
        "startPos": 3000000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 3800000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos75",
        "size": 5600000,
        "startPos": 6600000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 5100000,
        "startPos": 12200000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos75",
        "size": 1300000,
        "startPos": 17300000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 100000,
        "startPos": 18600000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos75",
        "size": 3900000,
        "startPos": 18700000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 22600000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos50",
        "size": 5000000,
        "startPos": 24600000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 29600000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos25",
        "size": 3100000,
        "startPos": 31300000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 3600000,
        "startPos": 34400000
    },
    {
        "chr": "10",
        "cytoBandType": "acen",
        "size": 2200000,
        "startPos": 38000000
    },
    {
        "chr": "10",
        "cytoBandType": "acen",
        "size": 2100000,
        "startPos": 40200000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 3800000,
        "startPos": 42300000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos25",
        "size": 3800000,
        "startPos": 46100000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 3000000,
        "startPos": 49900000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos100",
        "size": 8300000,
        "startPos": 52900000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 61200000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos100",
        "size": 6100000,
        "startPos": 64500000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 4300000,
        "startPos": 70600000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos50",
        "size": 2800000,
        "startPos": 74900000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 4300000,
        "startPos": 77700000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos100",
        "size": 5900000,
        "startPos": 82000000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 1600000,
        "startPos": 87900000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos75",
        "size": 3400000,
        "startPos": 89500000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 1200000,
        "startPos": 92900000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos50",
        "size": 2900000,
        "startPos": 94100000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 97000000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos50",
        "size": 2600000,
        "startPos": 99300000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 1100000,
        "startPos": 101900000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos25",
        "size": 1900000,
        "startPos": 103000000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 900000,
        "startPos": 104900000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos100",
        "size": 6100000,
        "startPos": 105800000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 3000000,
        "startPos": 111900000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos75",
        "size": 4200000,
        "startPos": 114900000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 119100000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos50",
        "size": 1400000,
        "startPos": 121700000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 123100000
    },
    {
        "chr": "10",
        "cytoBandType": "gpos50",
        "size": 3100000,
        "startPos": 127500000
    },
    {
        "chr": "10",
        "cytoBandType": "gneg",
        "size": 4934747,
        "startPos": 130600000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 0
    },
    {
        "chr": "11",
        "cytoBandType": "gpos50",
        "size": 7900000,
        "startPos": 2800000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 10700000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos50",
        "size": 3500000,
        "startPos": 12700000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 5500000,
        "startPos": 16200000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos100",
        "size": 4400000,
        "startPos": 21700000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 1100000,
        "startPos": 26100000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos75",
        "size": 3800000,
        "startPos": 27200000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 5400000,
        "startPos": 31000000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos100",
        "size": 7100000,
        "startPos": 36400000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 5300000,
        "startPos": 43500000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos75",
        "size": 2800000,
        "startPos": 48800000
    },
    {
        "chr": "11",
        "cytoBandType": "acen",
        "size": 2100000,
        "startPos": 51600000
    },
    {
        "chr": "11",
        "cytoBandType": "acen",
        "size": 2000000,
        "startPos": 53700000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos75",
        "size": 4200000,
        "startPos": 55700000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 59900000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos25",
        "size": 1700000,
        "startPos": 61700000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 2500000,
        "startPos": 63400000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos25",
        "size": 2500000,
        "startPos": 65900000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 68400000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos50",
        "size": 4800000,
        "startPos": 70400000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 75200000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos100",
        "size": 8500000,
        "startPos": 77100000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 85600000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos100",
        "size": 4500000,
        "startPos": 88300000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 92800000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos100",
        "size": 4900000,
        "startPos": 97200000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 800000,
        "startPos": 102100000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos100",
        "size": 7500000,
        "startPos": 102900000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 110400000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos50",
        "size": 2000000,
        "startPos": 112500000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 6700000,
        "startPos": 114500000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos50",
        "size": 2700000,
        "startPos": 121200000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 3900000,
        "startPos": 123900000
    },
    {
        "chr": "11",
        "cytoBandType": "gpos50",
        "size": 3000000,
        "startPos": 127800000
    },
    {
        "chr": "11",
        "cytoBandType": "gneg",
        "size": 4206516,
        "startPos": 130800000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 0
    },
    {
        "chr": "12",
        "cytoBandType": "gpos25",
        "size": 2100000,
        "startPos": 3300000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 4700000,
        "startPos": 5400000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos75",
        "size": 2700000,
        "startPos": 10100000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 12800000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos100",
        "size": 5200000,
        "startPos": 14800000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 1300000,
        "startPos": 20000000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos100",
        "size": 5200000,
        "startPos": 21300000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 1300000,
        "startPos": 26500000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos50",
        "size": 2900000,
        "startPos": 27800000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 30700000
    },
    {
        "chr": "12",
        "cytoBandType": "acen",
        "size": 2500000,
        "startPos": 33300000
    },
    {
        "chr": "12",
        "cytoBandType": "acen",
        "size": 2400000,
        "startPos": 35800000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos100",
        "size": 8200000,
        "startPos": 38200000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 46400000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos25",
        "size": 2400000,
        "startPos": 49100000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 51500000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos25",
        "size": 1700000,
        "startPos": 54900000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 1500000,
        "startPos": 56600000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos75",
        "size": 5000000,
        "startPos": 58100000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 63100000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos50",
        "size": 2600000,
        "startPos": 65100000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 3800000,
        "startPos": 67700000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos75",
        "size": 4200000,
        "startPos": 71500000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 4600000,
        "startPos": 75700000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos100",
        "size": 6400000,
        "startPos": 80300000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 86700000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos100",
        "size": 3600000,
        "startPos": 89000000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 3600000,
        "startPos": 92600000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos75",
        "size": 5400000,
        "startPos": 96200000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 101600000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos50",
        "size": 5200000,
        "startPos": 103800000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 109000000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos25",
        "size": 600000,
        "startPos": 111700000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 112300000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos50",
        "size": 2500000,
        "startPos": 114300000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 1300000,
        "startPos": 116800000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos50",
        "size": 2600000,
        "startPos": 118100000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 5200000,
        "startPos": 120700000
    },
    {
        "chr": "12",
        "cytoBandType": "gpos50",
        "size": 3400000,
        "startPos": 125900000
    },
    {
        "chr": "12",
        "cytoBandType": "gneg",
        "size": 4551895,
        "startPos": 129300000
    },
    {
        "chr": "13",
        "cytoBandType": "gvar",
        "size": 4500000,
        "startPos": 0
    },
    {
        "chr": "13",
        "cytoBandType": "stalk",
        "size": 5500000,
        "startPos": 4500000
    },
    {
        "chr": "13",
        "cytoBandType": "gvar",
        "size": 6300000,
        "startPos": 10000000
    },
    {
        "chr": "13",
        "cytoBandType": "acen",
        "size": 1600000,
        "startPos": 16300000
    },
    {
        "chr": "13",
        "cytoBandType": "acen",
        "size": 1600000,
        "startPos": 17900000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 3800000,
        "startPos": 19500000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos25",
        "size": 2200000,
        "startPos": 23300000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 25500000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos25",
        "size": 1100000,
        "startPos": 27800000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 28900000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos50",
        "size": 1800000,
        "startPos": 32200000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 1500000,
        "startPos": 34000000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos75",
        "size": 4600000,
        "startPos": 35500000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 5100000,
        "startPos": 40100000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos25",
        "size": 600000,
        "startPos": 45200000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 1500000,
        "startPos": 45800000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos50",
        "size": 3600000,
        "startPos": 47300000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 50900000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos100",
        "size": 4300000,
        "startPos": 55300000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 59600000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos75",
        "size": 3400000,
        "startPos": 62300000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 65700000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos100",
        "size": 4700000,
        "startPos": 68600000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 73300000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos50",
        "size": 1800000,
        "startPos": 75400000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 77200000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos100",
        "size": 8700000,
        "startPos": 79000000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 87700000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos100",
        "size": 5000000,
        "startPos": 90000000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 3200000,
        "startPos": 95000000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos25",
        "size": 1100000,
        "startPos": 98200000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 2400000,
        "startPos": 99300000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos100",
        "size": 3100000,
        "startPos": 101700000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 104800000
    },
    {
        "chr": "13",
        "cytoBandType": "gpos100",
        "size": 3300000,
        "startPos": 107000000
    },
    {
        "chr": "13",
        "cytoBandType": "gneg",
        "size": 4869878,
        "startPos": 110300000
    },
    {
        "chr": "14",
        "cytoBandType": "gvar",
        "size": 3700000,
        "startPos": 0
    },
    {
        "chr": "14",
        "cytoBandType": "stalk",
        "size": 4400000,
        "startPos": 3700000
    },
    {
        "chr": "14",
        "cytoBandType": "gvar",
        "size": 8000000,
        "startPos": 8100000
    },
    {
        "chr": "14",
        "cytoBandType": "acen",
        "size": 1500000,
        "startPos": 16100000
    },
    {
        "chr": "14",
        "cytoBandType": "acen",
        "size": 1500000,
        "startPos": 17600000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 5500000,
        "startPos": 19100000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos100",
        "size": 8700000,
        "startPos": 24600000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 33300000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos50",
        "size": 1300000,
        "startPos": 35300000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 1200000,
        "startPos": 36600000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos100",
        "size": 5700000,
        "startPos": 37800000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 3700000,
        "startPos": 43500000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos100",
        "size": 3700000,
        "startPos": 47200000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 3200000,
        "startPos": 50900000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos25",
        "size": 1400000,
        "startPos": 54100000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 55500000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos75",
        "size": 4000000,
        "startPos": 58100000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 62100000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos50",
        "size": 3100000,
        "startPos": 64800000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 67900000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos50",
        "size": 3600000,
        "startPos": 70200000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 5500000,
        "startPos": 73800000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos100",
        "size": 4300000,
        "startPos": 79300000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 1300000,
        "startPos": 83600000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos100",
        "size": 4900000,
        "startPos": 84900000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 89800000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos25",
        "size": 2800000,
        "startPos": 91900000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 1600000,
        "startPos": 94700000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos50",
        "size": 5100000,
        "startPos": 96300000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 101400000
    },
    {
        "chr": "14",
        "cytoBandType": "gpos50",
        "size": 800000,
        "startPos": 103200000
    },
    {
        "chr": "14",
        "cytoBandType": "gneg",
        "size": 3349540,
        "startPos": 104000000
    },
    {
        "chr": "15",
        "cytoBandType": "gvar",
        "size": 3900000,
        "startPos": 0
    },
    {
        "chr": "15",
        "cytoBandType": "stalk",
        "size": 4800000,
        "startPos": 3900000
    },
    {
        "chr": "15",
        "cytoBandType": "gvar",
        "size": 7100000,
        "startPos": 8700000
    },
    {
        "chr": "15",
        "cytoBandType": "acen",
        "size": 3200000,
        "startPos": 15800000
    },
    {
        "chr": "15",
        "cytoBandType": "acen",
        "size": 1700000,
        "startPos": 19000000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 5000000,
        "startPos": 20700000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos50",
        "size": 2400000,
        "startPos": 25700000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 28100000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos50",
        "size": 900000,
        "startPos": 30300000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 2400000,
        "startPos": 31200000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos75",
        "size": 6500000,
        "startPos": 33600000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 40100000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos25",
        "size": 800000,
        "startPos": 42800000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 1200000,
        "startPos": 43600000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos75",
        "size": 4700000,
        "startPos": 44800000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 49500000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos75",
        "size": 6200000,
        "startPos": 52900000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 200000,
        "startPos": 59100000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos25",
        "size": 4400000,
        "startPos": 59300000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 3500000,
        "startPos": 63700000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos25",
        "size": 100000,
        "startPos": 67200000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 200000,
        "startPos": 67300000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos25",
        "size": 5200000,
        "startPos": 67500000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 2500000,
        "startPos": 72700000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos25",
        "size": 1400000,
        "startPos": 75200000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 76600000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos50",
        "size": 3400000,
        "startPos": 78300000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 3500000,
        "startPos": 81700000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos50",
        "size": 3900000,
        "startPos": 85200000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 5200000,
        "startPos": 89100000
    },
    {
        "chr": "15",
        "cytoBandType": "gpos50",
        "size": 4200000,
        "startPos": 94300000
    },
    {
        "chr": "15",
        "cytoBandType": "gneg",
        "size": 4031392,
        "startPos": 98500000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 7900000,
        "startPos": 0
    },
    {
        "chr": "16",
        "cytoBandType": "gpos50",
        "size": 2600000,
        "startPos": 7900000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 10500000
    },
    {
        "chr": "16",
        "cytoBandType": "gpos50",
        "size": 2200000,
        "startPos": 12600000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 2000000,
        "startPos": 14800000
    },
    {
        "chr": "16",
        "cytoBandType": "gpos50",
        "size": 4400000,
        "startPos": 16800000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 3000000,
        "startPos": 21200000
    },
    {
        "chr": "16",
        "cytoBandType": "gpos50",
        "size": 3900000,
        "startPos": 24200000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 6500000,
        "startPos": 28100000
    },
    {
        "chr": "16",
        "cytoBandType": "acen",
        "size": 2000000,
        "startPos": 34600000
    },
    {
        "chr": "16",
        "cytoBandType": "acen",
        "size": 2000000,
        "startPos": 36600000
    },
    {
        "chr": "16",
        "cytoBandType": "gvar",
        "size": 8400000,
        "startPos": 38600000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 5600000,
        "startPos": 47000000
    },
    {
        "chr": "16",
        "cytoBandType": "gpos50",
        "size": 4100000,
        "startPos": 52600000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 700000,
        "startPos": 56700000
    },
    {
        "chr": "16",
        "cytoBandType": "gpos100",
        "size": 9300000,
        "startPos": 57400000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 4100000,
        "startPos": 66700000
    },
    {
        "chr": "16",
        "cytoBandType": "gpos50",
        "size": 2100000,
        "startPos": 70800000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 1200000,
        "startPos": 72900000
    },
    {
        "chr": "16",
        "cytoBandType": "gpos75",
        "size": 5100000,
        "startPos": 74100000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 2500000,
        "startPos": 79200000
    },
    {
        "chr": "16",
        "cytoBandType": "gpos50",
        "size": 2500000,
        "startPos": 81700000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 84200000
    },
    {
        "chr": "16",
        "cytoBandType": "gpos25",
        "size": 1600000,
        "startPos": 87100000
    },
    {
        "chr": "16",
        "cytoBandType": "gneg",
        "size": 1654753,
        "startPos": 88700000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 0
    },
    {
        "chr": "17",
        "cytoBandType": "gpos50",
        "size": 3200000,
        "startPos": 3300000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 4200000,
        "startPos": 6500000
    },
    {
        "chr": "17",
        "cytoBandType": "gpos75",
        "size": 5300000,
        "startPos": 10700000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 6200000,
        "startPos": 16000000
    },
    {
        "chr": "17",
        "cytoBandType": "acen",
        "size": 1800000,
        "startPos": 22200000
    },
    {
        "chr": "17",
        "cytoBandType": "acen",
        "size": 1800000,
        "startPos": 24000000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 6000000,
        "startPos": 25800000
    },
    {
        "chr": "17",
        "cytoBandType": "gpos50",
        "size": 6300000,
        "startPos": 31800000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 300000,
        "startPos": 38100000
    },
    {
        "chr": "17",
        "cytoBandType": "gpos25",
        "size": 2500000,
        "startPos": 38400000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 4000000,
        "startPos": 40900000
    },
    {
        "chr": "17",
        "cytoBandType": "gpos25",
        "size": 2500000,
        "startPos": 44900000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 47400000
    },
    {
        "chr": "17",
        "cytoBandType": "gpos75",
        "size": 7400000,
        "startPos": 50200000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 700000,
        "startPos": 57600000
    },
    {
        "chr": "17",
        "cytoBandType": "gpos75",
        "size": 2800000,
        "startPos": 58300000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 1500000,
        "startPos": 61100000
    },
    {
        "chr": "17",
        "cytoBandType": "gpos50",
        "size": 1600000,
        "startPos": 62600000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 64200000
    },
    {
        "chr": "17",
        "cytoBandType": "gpos75",
        "size": 3800000,
        "startPos": 67100000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 3900000,
        "startPos": 70900000
    },
    {
        "chr": "17",
        "cytoBandType": "gpos25",
        "size": 500000,
        "startPos": 74800000
    },
    {
        "chr": "17",
        "cytoBandType": "gneg",
        "size": 5895210,
        "startPos": 75300000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 0
    },
    {
        "chr": "18",
        "cytoBandType": "gpos50",
        "size": 4200000,
        "startPos": 2900000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 1400000,
        "startPos": 7100000
    },
    {
        "chr": "18",
        "cytoBandType": "gpos25",
        "size": 2400000,
        "startPos": 8500000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 4500000,
        "startPos": 10900000
    },
    {
        "chr": "18",
        "cytoBandType": "acen",
        "size": 1800000,
        "startPos": 15400000
    },
    {
        "chr": "18",
        "cytoBandType": "acen",
        "size": 1800000,
        "startPos": 17200000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 6000000,
        "startPos": 19000000
    },
    {
        "chr": "18",
        "cytoBandType": "gpos100",
        "size": 7700000,
        "startPos": 25000000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 4500000,
        "startPos": 32700000
    },
    {
        "chr": "18",
        "cytoBandType": "gpos75",
        "size": 6300000,
        "startPos": 37200000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 4700000,
        "startPos": 43500000
    },
    {
        "chr": "18",
        "cytoBandType": "gpos75",
        "size": 5600000,
        "startPos": 48200000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 2400000,
        "startPos": 53800000
    },
    {
        "chr": "18",
        "cytoBandType": "gpos50",
        "size": 2800000,
        "startPos": 56200000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 59000000
    },
    {
        "chr": "18",
        "cytoBandType": "gpos100",
        "size": 5200000,
        "startPos": 61600000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 66800000
    },
    {
        "chr": "18",
        "cytoBandType": "gpos25",
        "size": 4400000,
        "startPos": 68700000
    },
    {
        "chr": "18",
        "cytoBandType": "gneg",
        "size": 4977248,
        "startPos": 73100000
    },
    {
        "chr": "19",
        "cytoBandType": "gneg",
        "size": 6900000,
        "startPos": 0
    },
    {
        "chr": "19",
        "cytoBandType": "gpos25",
        "size": 7000000,
        "startPos": 6900000
    },
    {
        "chr": "19",
        "cytoBandType": "gneg",
        "size": 100000,
        "startPos": 13900000
    },
    {
        "chr": "19",
        "cytoBandType": "gpos25",
        "size": 2300000,
        "startPos": 14000000
    },
    {
        "chr": "19",
        "cytoBandType": "gneg",
        "size": 3700000,
        "startPos": 16300000
    },
    {
        "chr": "19",
        "cytoBandType": "gvar",
        "size": 4400000,
        "startPos": 20000000
    },
    {
        "chr": "19",
        "cytoBandType": "acen",
        "size": 2100000,
        "startPos": 24400000
    },
    {
        "chr": "19",
        "cytoBandType": "acen",
        "size": 2100000,
        "startPos": 26500000
    },
    {
        "chr": "19",
        "cytoBandType": "gvar",
        "size": 3800000,
        "startPos": 28600000
    },
    {
        "chr": "19",
        "cytoBandType": "gneg",
        "size": 3100000,
        "startPos": 32400000
    },
    {
        "chr": "19",
        "cytoBandType": "gpos25",
        "size": 2800000,
        "startPos": 35500000
    },
    {
        "chr": "19",
        "cytoBandType": "gneg",
        "size": 400000,
        "startPos": 38300000
    },
    {
        "chr": "19",
        "cytoBandType": "gpos25",
        "size": 4700000,
        "startPos": 38700000
    },
    {
        "chr": "19",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 43400000
    },
    {
        "chr": "19",
        "cytoBandType": "gpos25",
        "size": 2800000,
        "startPos": 45200000
    },
    {
        "chr": "19",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 48000000
    },
    {
        "chr": "19",
        "cytoBandType": "gpos25",
        "size": 2200000,
        "startPos": 51400000
    },
    {
        "chr": "19",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 53600000
    },
    {
        "chr": "19",
        "cytoBandType": "gpos25",
        "size": 2828983,
        "startPos": 56300000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 5100000,
        "startPos": 0
    },
    {
        "chr": "20",
        "cytoBandType": "gpos75",
        "size": 4100000,
        "startPos": 5100000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 2900000,
        "startPos": 9200000
    },
    {
        "chr": "20",
        "cytoBandType": "gpos75",
        "size": 5800000,
        "startPos": 12100000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 17900000
    },
    {
        "chr": "20",
        "cytoBandType": "gpos25",
        "size": 1000000,
        "startPos": 21300000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 22300000
    },
    {
        "chr": "20",
        "cytoBandType": "acen",
        "size": 1900000,
        "startPos": 25600000
    },
    {
        "chr": "20",
        "cytoBandType": "acen",
        "size": 1900000,
        "startPos": 27500000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 2700000,
        "startPos": 29400000
    },
    {
        "chr": "20",
        "cytoBandType": "gpos25",
        "size": 2300000,
        "startPos": 32100000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 3200000,
        "startPos": 34400000
    },
    {
        "chr": "20",
        "cytoBandType": "gpos75",
        "size": 4100000,
        "startPos": 37600000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 400000,
        "startPos": 41700000
    },
    {
        "chr": "20",
        "cytoBandType": "gpos25",
        "size": 4300000,
        "startPos": 42100000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 46400000
    },
    {
        "chr": "20",
        "cytoBandType": "gpos75",
        "size": 5200000,
        "startPos": 49800000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 1500000,
        "startPos": 55000000
    },
    {
        "chr": "20",
        "cytoBandType": "gpos50",
        "size": 1900000,
        "startPos": 56500000
    },
    {
        "chr": "20",
        "cytoBandType": "gneg",
        "size": 4625520,
        "startPos": 58400000
    },
    {
        "chr": "21",
        "cytoBandType": "gvar",
        "size": 2800000,
        "startPos": 0
    },
    {
        "chr": "21",
        "cytoBandType": "stalk",
        "size": 4000000,
        "startPos": 2800000
    },
    {
        "chr": "21",
        "cytoBandType": "gvar",
        "size": 4100000,
        "startPos": 6800000
    },
    {
        "chr": "21",
        "cytoBandType": "acen",
        "size": 2300000,
        "startPos": 10900000
    },
    {
        "chr": "21",
        "cytoBandType": "acen",
        "size": 1100000,
        "startPos": 13200000
    },
    {
        "chr": "21",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 14300000
    },
    {
        "chr": "21",
        "cytoBandType": "gpos100",
        "size": 7600000,
        "startPos": 16400000
    },
    {
        "chr": "21",
        "cytoBandType": "gneg",
        "size": 2800000,
        "startPos": 24000000
    },
    {
        "chr": "21",
        "cytoBandType": "gpos75",
        "size": 4700000,
        "startPos": 26800000
    },
    {
        "chr": "21",
        "cytoBandType": "gneg",
        "size": 4300000,
        "startPos": 31500000
    },
    {
        "chr": "21",
        "cytoBandType": "gpos50",
        "size": 2000000,
        "startPos": 35800000
    },
    {
        "chr": "21",
        "cytoBandType": "gneg",
        "size": 1900000,
        "startPos": 37800000
    },
    {
        "chr": "21",
        "cytoBandType": "gpos50",
        "size": 2900000,
        "startPos": 39700000
    },
    {
        "chr": "21",
        "cytoBandType": "gneg",
        "size": 5529895,
        "startPos": 42600000
    },
    {
        "chr": "22",
        "cytoBandType": "gvar",
        "size": 3800000,
        "startPos": 0
    },
    {
        "chr": "22",
        "cytoBandType": "stalk",
        "size": 4500000,
        "startPos": 3800000
    },
    {
        "chr": "22",
        "cytoBandType": "gvar",
        "size": 3900000,
        "startPos": 8300000
    },
    {
        "chr": "22",
        "cytoBandType": "acen",
        "size": 2500000,
        "startPos": 12200000
    },
    {
        "chr": "22",
        "cytoBandType": "acen",
        "size": 3200000,
        "startPos": 14700000
    },
    {
        "chr": "22",
        "cytoBandType": "gneg",
        "size": 4300000,
        "startPos": 17900000
    },
    {
        "chr": "22",
        "cytoBandType": "gpos25",
        "size": 1300000,
        "startPos": 22200000
    },
    {
        "chr": "22",
        "cytoBandType": "gneg",
        "size": 2400000,
        "startPos": 23500000
    },
    {
        "chr": "22",
        "cytoBandType": "gpos50",
        "size": 3700000,
        "startPos": 25900000
    },
    {
        "chr": "22",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 29600000
    },
    {
        "chr": "22",
        "cytoBandType": "gpos50",
        "size": 5400000,
        "startPos": 32200000
    },
    {
        "chr": "22",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 37600000
    },
    {
        "chr": "22",
        "cytoBandType": "gpos50",
        "size": 3200000,
        "startPos": 41000000
    },
    {
        "chr": "22",
        "cytoBandType": "gneg",
        "size": 4200000,
        "startPos": 44200000
    },
    {
        "chr": "22",
        "cytoBandType": "gpos50",
        "size": 1000000,
        "startPos": 48400000
    },
    {
        "chr": "22",
        "cytoBandType": "gneg",
        "size": 1904566,
        "startPos": 49400000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 4300000,
        "startPos": 0
    },
    {
        "chr": "X",
        "cytoBandType": "gpos50",
        "size": 1700000,
        "startPos": 4300000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 3500000,
        "startPos": 6000000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos50",
        "size": 7600000,
        "startPos": 9500000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 17100000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos50",
        "size": 2600000,
        "startPos": 19300000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 3000000,
        "startPos": 21900000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos100",
        "size": 4400000,
        "startPos": 24900000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 2200000,
        "startPos": 29300000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos100",
        "size": 6100000,
        "startPos": 31500000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 4800000,
        "startPos": 37600000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos75",
        "size": 4000000,
        "startPos": 42400000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 3400000,
        "startPos": 46400000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos25",
        "size": 5000000,
        "startPos": 49800000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 3300000,
        "startPos": 54800000
    },
    {
        "chr": "X",
        "cytoBandType": "acen",
        "size": 2500000,
        "startPos": 58100000
    },
    {
        "chr": "X",
        "cytoBandType": "acen",
        "size": 2400000,
        "startPos": 60600000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 1600000,
        "startPos": 63000000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos50",
        "size": 3200000,
        "startPos": 64600000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 4000000,
        "startPos": 67800000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos50",
        "size": 2100000,
        "startPos": 71800000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 2100000,
        "startPos": 73900000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos100",
        "size": 8600000,
        "startPos": 76000000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 1600000,
        "startPos": 84600000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos100",
        "size": 5600000,
        "startPos": 86200000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 91800000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos75",
        "size": 4800000,
        "startPos": 93500000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 4300000,
        "startPos": 98300000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos50",
        "size": 1100000,
        "startPos": 102600000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 5000000,
        "startPos": 103700000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos75",
        "size": 7800000,
        "startPos": 108700000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 116500000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos100",
        "size": 7800000,
        "startPos": 120900000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 128700000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos25",
        "size": 3200000,
        "startPos": 130400000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 4400000,
        "startPos": 133600000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos75",
        "size": 2300000,
        "startPos": 138000000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 1800000,
        "startPos": 140300000
    },
    {
        "chr": "X",
        "cytoBandType": "gpos100",
        "size": 5000000,
        "startPos": 142100000
    },
    {
        "chr": "X",
        "cytoBandType": "gneg",
        "size": 8170560,
        "startPos": 147100000
    },
    {
        "chr": "Y",
        "cytoBandType": "gneg",
        "size": 2500000,
        "startPos": 0
    },
    {
        "chr": "Y",
        "cytoBandType": "gpos50",
        "size": 500000,
        "startPos": 2500000
    },
    {
        "chr": "Y",
        "cytoBandType": "gneg",
        "size": 8600000,
        "startPos": 3000000
    },
    {
        "chr": "Y",
        "cytoBandType": "acen",
        "size": 900000,
        "startPos": 11600000
    },
    {
        "chr": "Y",
        "cytoBandType": "acen",
        "size": 900000,
        "startPos": 12500000
    },
    {
        "chr": "Y",
        "cytoBandType": "gneg",
        "size": 1700000,
        "startPos": 13400000
    },
    {
        "chr": "Y",
        "cytoBandType": "gpos50",
        "size": 4700000,
        "startPos": 15100000
    },
    {
        "chr": "Y",
        "cytoBandType": "gneg",
        "size": 2300000,
        "startPos": 19800000
    },
    {
        "chr": "Y",
        "cytoBandType": "gpos50",
        "size": 4100000,
        "startPos": 22100000
    },
    {
        "chr": "Y",
        "cytoBandType": "gneg",
        "size": 2600000,
        "startPos": 26200000
    },
    {
        "chr": "Y",
        "cytoBandType": "gvar",
        "size": 30573566,
        "startPos": 28800000
    }
];
