use itertools::Itertools;
use std::{
    collections::{HashMap, HashSet},
    hash::Hash,
};
use ndarray;

/// Flatten site array, i.e. treat individuals as haploid
pub fn flatten_site(site: Vec<Vec<u32>>) -> Vec<u32> {
    site.into_iter().flatten().collect()
}

/// Get SFS entry (i,j) from single site of genotypes (biallelic only)
pub fn site_to_entry(site: Vec<u32>, sample_map: &HashMap<usize, String>) -> Vec<usize> {
    let mut populations: Vec<&String> = Vec::from_iter(sample_map.values().sorted());
    populations.dedup(); // make unique list of populations
    let mut ntons: Vec<usize> = vec![];

    for population in populations {
        let nton: &usize = &site
            .iter()
            .enumerate()
            .filter(|(idx, _)| *sample_map.get(idx).unwrap() == *population)
            .filter(|(_, val)| **val != 0 as u32)
            .count();

        ntons.push(*nton)
    }

    ntons
}

/// Compute bSFS from calls; return indices of entries in bSFS matrix
pub fn bsfs_indices(
    calls: Vec<Vec<Vec<u32>>>,
    sample_map: HashMap<usize, String>,
) -> Vec<Vec<usize>> {
    let indices = calls
        .into_iter()
        .map(|site| flatten_site(site))
        .map(|site| site_to_entry(site, &sample_map))
        .collect();

    indices
}

/// Count number of haplotypes per population
pub fn n_haps_per_pop(sample_map: &HashMap<usize, String>) -> HashMap<String, usize> {
    let mut populations: Vec<&String> = Vec::from_iter(sample_map.values().sorted());
    populations.dedup(); // make unique list of populations

    let mut counts: HashMap<String, usize> = Default::default();

    for individual in sample_map.keys() {
        let pop: String = if let Some(popname) = sample_map.get(individual) {
            popname.to_string()
        } else {
            continue;
        };

        *counts.entry(pop).or_default() += 1;
    }

    counts
}

/// Get bSFS matrix for block
pub fn bsfs_matrix(calls: Vec<Vec<Vec<u32>>>, sample_map: HashMap<usize, String>) -> Vec<Vec<u32>> {
    // Count number of haplotypes
    let mut haps_per_pop = n_haps_per_pop(&sample_map)

    let bsfs_matrix: ndarray::ArrayBase<> = ndarray::ArrayBase::zeros((&haps_per_pop[0], &haps_per_pop[1]));
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use rstest::{fixture, rstest};

    #[fixture]
    fn single_call_arr() -> Vec<Vec<Vec<u32>>> {
        let single_call_arr: Vec<Vec<Vec<u32>>> =
            vec![vec![vec![0, 0], vec![0, 1], vec![0, 0], vec![0, 1]]];
        single_call_arr
    }

    #[fixture]
    fn flattened_single_call() -> Vec<u32> {
        vec![0, 0, 0, 1, 0, 0, 0, 1]
    }

    #[fixture]
    fn block_calls() -> Vec<Vec<Vec<u32>>> {
        vec![
            vec![vec![0, 1], vec![1, 0], vec![0, 0], vec![1, 1]],
            vec![vec![1, 1], vec![1, 1], vec![1, 1], vec![1, 1]],
            vec![vec![0, 0], vec![0, 1], vec![0, 0], vec![0, 0]],
        ]
    }

    #[fixture]
    fn sample_map() -> HashMap<usize, String> {
        let mut samplemap: HashMap<usize, String> = HashMap::new();
        samplemap.insert(0, "popA".to_string());
        samplemap.insert(1, "popA".to_string());
        samplemap.insert(2, "popA".to_string());
        samplemap.insert(3, "popA".to_string());
        samplemap.insert(4, "popB".to_string());
        samplemap.insert(5, "popB".to_string());
        samplemap.insert(6, "popB".to_string());
        samplemap.insert(7, "popB".to_string());

        samplemap
    }

    #[rstest]
    fn test_flatten_site(single_call_arr: Vec<Vec<Vec<u32>>>, flattened_single_call: Vec<u32>) {
        let flattened: Vec<u32> = flatten_site(single_call_arr[0].clone());
        let expected: Vec<u32> = flattened_single_call;

        assert_eq!(flattened, expected)
    }

    #[rstest]
    fn test_site_to_entry(flattened_single_call: Vec<u32>, sample_map: HashMap<usize, String>) {
        let entry: Vec<usize> = site_to_entry(flattened_single_call, &sample_map);
        let expected = vec![1, 1];

        assert_eq!(entry, expected)
    }

    #[rstest]
    fn test_bsfs_indices(block_calls: Vec<Vec<Vec<u32>>>, sample_map: HashMap<usize, String>) {
        let bsfs_indices = bsfs_indices(block_calls, sample_map);
        let expected = vec![[2, 2], [4, 4], [1, 0]];

        assert_eq!(bsfs_indices, expected)
    }

    #[rstest]
    fn test_bsfs_matrix(block_calls: Vec<Vec<Vec<u32>>>, sample_map: HashMap<usize, String>) {
        let block_bsfs = bsfs_matrix(block_calls, sample_map);
        let expected = vec![[0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];

        assert_eq!(block_bsfs, expected)
    }

    #[rstest]
    fn test_n_haps_per_pop(sample_map: HashMap<usize, String>) {
        let expected = HashMap::from([("popA".to_string(), 4), ("popB".to_string(), 4)]);

        let counts = n_haps_per_pop(&sample_map);

        assert_eq!(counts, expected)
    }
}
