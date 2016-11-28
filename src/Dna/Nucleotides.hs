{-# LANGUAGE FlexibleInstances, RankNTypes #-}
module Dna.Nucleotides (
    Nucleotide, Dna, DnaList, Rna, (+:+), isDnaPrefixOf, lengthDna,
    groupDna, tailDna, readsDna, DnaGraph, nullDna, liftDna2, liftDna, takeDna,
    overlapGraphFromNode, initOverlaGraph, dnaToRna, dropDna, dnaComplement)
where

import Data.List
import Data.Functor
import Data.STRef
import Control.Monad
import Control.Applicative
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data
import Numeric.LinearAlgebra.Devel
import Control.Monad.ST

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


newtype Dna = DnaIntern [Nucleotide]
    deriving (Ord, Eq)
newtype Rna = RnaIntern [Nucleotide]
    deriving (Ord, Eq)
type DnaList = [Dna]

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


data DnaGraph = DnaGraphInternal {
    node :: S.Set Dna,
    connection :: Matrix I
}

isConnected :: Dna -> Dna -> Bool
isConnected (DnaIntern first) (DnaIntern second) = and $ zipWith (==) (tail first) second

initOverlaGraph :: DnaList -> (forall s. STMatrix s I -> (Dna -> Int) -> ST s ()) -> DnaGraph
initOverlaGraph nodes f = DnaGraphInternal setNode (runSTMatrix $ do
        m <- newMatrix 0 lenNodes lenNodes
        f m (`S.findIndex` setNode)
        return m
    )
    where
        setNode = S.fromList nodes
        lenNodes  = length setNode

overlapGraphFromNode nodes = DnaGraphInternal seqNode
            (build (lenMat,lenMat) (\first second ->
                if isConnected (S.elemAt (fromIntegral first) seqNode ) (S.elemAt (fromIntegral second) seqNode )
                    then 1
                    else 0))
        where seqNode = S.fromList nodes
              lenMat  = length seqNode

instance Show DnaGraph where
    show graph = intercalate "\n"
            [link | row <- [0 .. rows (connection graph) - 1],
                    any (isLink row) [0 .. cols (connection graph) - 1],
                    link <- [show (S.elemAt row (node graph)) ++ " -> " ++
                             intercalate "," (concatMap (\col -> replicate (fromIntegral $ linkCount row col) $ show  $ S.elemAt col (node graph))
                                                (filter (isLink row) [0 .. cols (connection graph) - 1]))]]
                where
                    linkCount row col = connection graph ! row ! col
                    isLink row col =   linkCount row col >= 1

mapNucleotidDnaRna T = U
mapNucleotidDnaRna a = a

dnaToRna :: Dna -> Rna
dnaToRna = RnaIntern . map mapNucleotidDnaRna . unDna
