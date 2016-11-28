module Dna.MetaInfo (
    nucleotideFrequency
)
where

import Dna.Nucleotides
import Data.List

nucleotideFrequency :: Dna -> [Int]
nucleotideFrequency dna = map lengthDna $ groupDna $ liftDna sort dna
