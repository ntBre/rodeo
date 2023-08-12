use std::ops::{Index, IndexMut};

#[derive(Clone, Copy)]
pub(crate) struct VertexDescriptor(pub(crate) usize);

#[derive(Clone, Copy)]
pub(crate) struct EdgeDescriptor(pub(crate) usize);

#[derive(Default)]
pub(crate) struct Graph<V, E> {
    pub(crate) vertices: Vec<V>,
    pub(crate) edges: Vec<E>,
}

impl<V, E> Graph<V, E>
where
    V: Default,
    E: Default,
{
    /// allocate a new vertex and return its index
    pub(crate) fn add_vertex(&mut self) -> VertexDescriptor {
        let next = self.vertices.len();
        self.vertices.push(V::default());
        VertexDescriptor(next)
    }

    /// allocate a new edge and return its index
    pub(crate) fn add_edge(&mut self) -> EdgeDescriptor {
        let next = self.edges.len();
        self.edges.push(E::default());
        EdgeDescriptor(next)
    }
}

impl<V, E> Index<VertexDescriptor> for Graph<V, E> {
    type Output = V;

    fn index(&self, index: VertexDescriptor) -> &Self::Output {
        &self.vertices[index.0]
    }
}

impl<V, E> IndexMut<VertexDescriptor> for Graph<V, E> {
    fn index_mut(&mut self, index: VertexDescriptor) -> &mut Self::Output {
        &mut self.vertices[index.0]
    }
}

impl<V, E> Index<EdgeDescriptor> for Graph<V, E> {
    type Output = E;

    fn index(&self, index: EdgeDescriptor) -> &Self::Output {
        &self.edges[index.0]
    }
}

impl<V, E> IndexMut<EdgeDescriptor> for Graph<V, E> {
    fn index_mut(&mut self, index: EdgeDescriptor) -> &mut Self::Output {
        &mut self.edges[index.0]
    }
}
