{-# LANGUAGE FlexibleInstances, RankNTypes #-}
module Dna.Nucleotides (
    Nucleotide, Dna, DnaList, Rna, (+:+),
    isDnaPrefixOf, lengthDna,
    groupDna, tailDna, readsDna, nullDna,
    liftDna2, liftDna, takeDna, dnaToRna,
    dropDna, dnaComplement, FastaSequence,
    parseFastaFile, highestGC)
where

import           Text.Printf
import           Data.Ratio
import           Control.Monad.IO.Class
import           Data.List
import           Data.Functor
import           Data.STRef
import           Control.Monad
import           Control.Applicative
import           Numeric.LinearAlgebra
import           Numeric.LinearAlgebra.Data
import           Numeric.LinearAlgebra.Devel
import           Control.Monad.ST
import qualified Text.Parser.Combinators as P
import qualified Text.Parser.Char as P
import qualified Text.Parser.Token as P
import qualified Text.Trifecta.Parser as P


import qualified Data.ByteString.Lazy.Char8 as C
import qualified Data.Set as S


data Nucleotide = A | C | G | T | U
    deriving (Eq, Ord)

instance Show Nucleotide where
    show A = "A"
    show C = "C"
    show G = "G"
    show T = "T"
    show U = "U"

data FastaSequence a = Fasta {
    _fastaId :: String,
    _fastaInfo :: a}

class FastaParseInfo a where
    parse :: P.TokenParsing m => m a

instance Functor FastaSequence where
    fmap f (Fasta ident datum) = Fasta ident (f datum)

instance FastaParseInfo Dna where
    parse = DnaIntern . concat <$> P.sepEndBy (many parseDnaNucleotide) P.newline

instance FastaParseInfo Rna where
    parse = RnaIntern . concat <$> P.sepEndBy (many parseRnaNucleotide) P.newline

newtype Dna = DnaIntern [Nucleotide]
    deriving (Ord, Eq)
newtype Rna = RnaIntern [Nucleotide]
    deriving (Ord, Eq)
type DnaList = [Dna]

parseFastaFile :: (FastaParseInfo a, MonadIO m) => String -> m (Maybe [FastaSequence a])
parseFastaFile = P.parseFromFile parseFastaList

parseDnaNucleotide :: P.TokenParsing m => m Nucleotide
parseDnaNucleotide = P.choice [P.string "A" *> pure A, P.string "T" *> pure T, P.string "G" *> pure G, P.string "C" *> pure C]

parseRnaNucleotide :: P.TokenParsing m => m Nucleotide
parseRnaNucleotide = P.choice [P.string "A" *> pure A, P.string "U" *> pure U, P.string "G" *> pure G, P.string "C" *> pure C]


parseFastaList :: (FastaParseInfo a, P.TokenParsing m) => m [FastaSequence a]
parseFastaList = P.many parseFasta

parseFasta :: (FastaParseInfo a, P.TokenParsing m) => m (FastaSequence a)
parseFasta = Fasta
                <$> (P.string ">" *> P.manyTill P.anyChar P.newline)
                <*> parse

unDna (DnaIntern dna) = dna
unRna (DnaIntern rna) = rna

liftDna :: ([Nucleotide] -> [Nucleotide]) -> Dna -> Dna
liftDna f (DnaIntern dna) = DnaIntern (f dna)

liftDna2 :: ([Nucleotide] -> [Nucleotide] -> a) -> Dna -> Dna -> a
liftDna2 f (DnaIntern dna1) (DnaIntern dna2) =  f dna1 dna2

groupDna :: Dna -> DnaList
groupDna = map DnaIntern . group . unDna

dropDna :: Int -> Dna -> Dna
dropDna n = liftDna $ drop n

takeDna :: Int -> Dna -> Dna
takeDna n = liftDna $ take n

tailDna :: Dna -> Dna
tailDna = liftDna tail

lengthDna :: Dna -> Int
lengthDna (DnaIntern dna) = length dna

isDnaPrefixOf :: Dna -> Dna -> Bool
isDnaPrefixOf  = liftDna2 isPrefixOf

(+:+) :: Dna -> Dna -> Dna
l +:+ r = DnaIntern $ liftDna2 (++) l r

nullDna :: Dna
nullDna = DnaIntern []

instance Show Dna where
    show (DnaIntern dna) = foldr ((++) . show) "" dna

instance Show Rna where
    show (RnaIntern rna) = foldr ((++) . show) "" rna

instance (Show a) => Show (FastaSequence a) where
    show (Fasta ident datum) = ">" ++ ident ++ "\n" ++ show datum ++ "\n"

instance {-# OVERLAPPING #-} Show (FastaSequence Rational) where
    show (Fasta ident datum) = ">" ++ ident ++ "\n" ++ show (fromRational datum :: Double)  ++ "\n"

instance {-# OVERLAPPING #-} Show DnaList where
    show = intercalate "\n" . map show

readsNucleotide          :: Int -> ReadS Nucleotide
readsNucleotide _ ('A':xs) = [(A,xs)]
readsNucleotide _ ('C':xs) = [(C,xs)]
readsNucleotide _ ('G':xs) = [(G,xs)]
readsNucleotide _ ('T':xs) = [(T,xs)]
readsNucleotide _ ('U':xs) = [(U,xs)]


instance Read Nucleotide where
    readsPrec = readsNucleotide

highestGC :: [FastaSequence Dna] -> FastaSequence Rational
highestGC = last . sortOn  _fastaInfo . fmap  ((`frequency` [C,G]) <$>)

readsDna        :: Int -> ReadS Dna
readsDna  _ []  = [(DnaIntern [],[])]
readsDna  d s@(x:xs)
    | x `elem` ['A','C','G','T'] = [(DnaIntern $ read [x]:dna, rest) | (DnaIntern dna,rest) <- readsDna (d+1) xs ]
    | otherwise = [(DnaIntern [],s)]

readsRna        :: Int -> ReadS Rna
readsRna  _ []  = [(RnaIntern [],[])]
readsRna  d s@(x:xs)
    | x `elem` ['A','C','G','U'] = [(RnaIntern $ read [x]:rna, rest) | (RnaIntern rna,rest) <- readsRna (d+1) xs ]
    | otherwise = [(RnaIntern [],s)]

readsDnaList    :: Int -> ReadS DnaList
readsDnaList _ [] = [([],[])]
readsDnaList d ('\n':xs) = readsDnaList (d+1) xs
readsDnaList d s@(x:xs)
    | x `notElem` ['A','C','G','T'] = [([],s)]
    | otherwise = [(DnaIntern dna:restDna,rest) |
            (DnaIntern dna,restRead) <- readsDna d s,
            (restDna,rest) <- readsDnaList (d+length dna) restRead]

mapComplementDan A = T
mapComplementDan T = A
mapComplementDan C = G
mapComplementDan G = C

dnaComplement :: Dna -> Dna
dnaComplement = liftDna $ reverse . map mapComplementDan

instance Read Dna where
    readsPrec = readsDna

instance Read Rna where
    readsPrec = readsRna

instance {-# OVERLAPPING #-} Read DnaList where
    readsPrec = readsDnaList


isConnected :: Dna -> Dna -> Bool
isConnected (DnaIntern first) (DnaIntern second) = and $ zipWith (==) (tail first) second

dnaToRna :: Dna -> Rna
dnaToRna = RnaIntern . map mapNucleotidDnaRna . unDna
    where mapNucleotidDnaRna T = U
          mapNucleotidDnaRna a = a

frequency :: Dna -> [Nucleotide] -> Rational
frequency (DnaIntern d) pop = fromIntegral (length (filter (`elem` pop) d)) %  fromIntegral (length d)
