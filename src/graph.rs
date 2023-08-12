use std::ops::{Index, IndexMut};

#[derive(Default)]
pub(crate) struct Graph<T> {
    vertices: Vec<T>,
}

impl<T> Graph<T>
where
    T: Default,
{
    /// allocate a new vertex and return its index
    pub(crate) fn add_vertex(&mut self) -> usize {
        let next = self.vertices.len();
        self.vertices.push(T::default());
        next
    }

    /// insert a new edge between
    pub(crate) fn add_edge(&self, atom1: usize, atom2: usize) -> usize {
        todo!()
    }
}

impl<T> Index<usize> for Graph<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.vertices[index]
    }
}

impl<T> IndexMut<usize> for Graph<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.vertices[index]
    }
}
