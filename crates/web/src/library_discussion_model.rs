//! One exact revision's human discussion state, independent of browser transport.
use collagene_library_protocol::{self as wire, discussion as api};
#[derive(Debug, Clone)]
pub struct Review {
    pub target: api::Target,
    pub generation: u64,
    pub events: Vec<api::Event>,
    pub member: api::Member,
    pub next_offset: Option<u32>,
    pub fresh: bool,
}
impl Review {
    pub fn page(
        previous: Option<&Self>,
        request: &api::RequestEnvelope,
        response: &api::ResponseEnvelope,
        expected_lab: &wire::LabId,
    ) -> Result<Self, String> {
        response
            .validate_for(request, expected_lab)
            .map_err(|e| e.to_string())?;
        let api::Request::List { target, offset, .. } = &request.request else {
            return Err("Expected a discussion listing.".into());
        };
        let api::Outcome::Success(api::Success::Listed(reply)) = &response.outcome else {
            return Err("Expected a discussion page.".into());
        };
        let offset = *offset;
        let member = response.member.clone();
        let mut events = if offset == 0 {
            Vec::new()
        } else {
            let previous = previous.ok_or("Reload discussion first.")?;
            if previous.target != *target
                || previous.events.len() != offset as usize
                || previous.member != member
                || previous.generation != reply.generation
                || previous.next_offset != Some(offset)
                || !previous.fresh
            {
                return Err("Discussion context changed; reload it.".into());
            }
            previous.events.clone()
        };
        events.extend(reply.events.as_slice().iter().cloned());
        Ok(Self {
            target: target.clone(),
            generation: reply.generation,
            events,
            member,
            next_offset: reply.next_offset,
            fresh: true,
        })
    }
    pub fn can_write(&self) -> bool {
        self.fresh
            && self.next_offset.is_none()
            && self.member.role != api::Role::Reader
            && self.generation < api::MAX_EVENTS
    }
    pub fn append_request(&self, action: api::Action) -> Result<api::Request, String> {
        if !self.can_write() {
            return Err("Load the complete current discussion before posting.".into());
        }
        action.validate().map_err(|e| e.to_string())?;
        Ok(api::Request::Append {
            target: self.target.clone(),
            expected_generation: self.generation,
            expected_actor: self.member.author.id.clone(),
            action,
        })
    }
    pub fn resolved(&self, thread: u64) -> Option<bool> {
        let root = self.events.iter().find(|e| e.sequence == thread)?;
        if !matches!(root.action, api::Action::OpenThread { .. }) {
            return None;
        }
        let mut resolved = false;
        for e in &self.events {
            if let api::Action::SetResolved {
                thread: t,
                resolved: r,
            } = e.action
            {
                if t == thread {
                    resolved = r;
                }
            }
        }
        Some(resolved)
    }
    pub fn can_resolve(&self, thread: u64) -> bool {
        self.can_write()
            && self.events.iter().any(|e| {
                e.sequence == thread
                    && matches!(e.action, api::Action::OpenThread { .. })
                    && (self.member.role == api::Role::Maintainer
                        || e.author.id == self.member.author.id)
            })
    }
    pub fn can_reply(&self, thread: u64) -> bool {
        self.can_write() && self.resolved(thread) == Some(false)
    }
    pub fn latest_reviews(&self) -> Vec<api::Event> {
        let mut latest = std::collections::BTreeMap::new();
        for event in &self.events {
            if matches!(event.action, api::Action::Review { .. }) {
                latest.insert(event.author.id.as_str(), event.clone());
            }
        }
        latest.into_values().collect()
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (api::Target, api::Member) {
        (
            api::Target {
                document: wire::DocumentId::new(format!("doc-{}", "a".repeat(32))).unwrap(),
                revision: wire::RevisionId::new("b".repeat(64)).unwrap(),
            },
            api::Member {
                author: api::Author {
                    id: wire::Text::new("scientist").unwrap(),
                    display_name: wire::Text::new("Scientist").unwrap(),
                },
                role: api::Role::Contributor,
            },
        )
    }
    fn event(n: u64, member: &api::Member, action: api::Action) -> api::Event {
        api::Event {
            sequence: n,
            author: member.author.clone(),
            role: member.role,
            created_at: 1,
            action,
        }
    }
    fn page(
        previous: Option<&Review>,
        target: &api::Target,
        offset: u32,
        page: api::Page,
        member: api::Member,
    ) -> Result<Review, String> {
        let lab = wire::LabId::new(format!("lab-{}", "0".repeat(32))).unwrap();
        let request = api::RequestEnvelope {
            request_id: 9,
            request: api::Request::List {
                target: target.clone(),
                offset,
                limit: 200,
                expected_generation: previous.map(|p| p.generation),
            },
        };
        let response = api::ResponseEnvelope {
            request_id: 9,
            context: wire::LabContext {
                lab_id: lab.clone(),
                name: wire::Text::new("Lab").unwrap(),
            },
            member,
            outcome: api::Outcome::Success(api::Success::Listed(page)),
        };
        Review::page(previous, &request, &response, &lab)
    }
    #[test]
    fn complete_snapshot_and_fresh_revision_required_before_posting() {
        let (target, member) = fixture();
        let events = (1..=200)
            .map(|n| {
                event(
                    n,
                    &member,
                    api::Action::OpenThread {
                        body: wire::Text::new("thread").unwrap(),
                    },
                )
            })
            .collect();
        let first = page(
            None,
            &target,
            0,
            api::Page {
                target: target.clone(),
                generation: 201,
                events: wire::Rows::new(events).unwrap(),
                next_offset: Some(200),
            },
            member.clone(),
        )
        .unwrap();
        assert!(!first.can_write());
        let last = event(
            201,
            &member,
            api::Action::Review {
                decision: api::ReviewDecision::Approved,
                body: wire::Text::new("").unwrap(),
            },
        );
        let last_page = api::Page {
            target: target.clone(),
            generation: 201,
            events: wire::Rows::new(vec![last]).unwrap(),
            next_offset: None,
        };
        let mut full = page(
            Some(&first),
            &target,
            200,
            last_page.clone(),
            member.clone(),
        )
        .unwrap();
        assert!(full.can_write());
        assert_eq!(full.latest_reviews().len(), 1);
        full.fresh = false;
        assert!(!full.can_write());
        assert!(
            full.append_request(api::Action::Review {
                decision: api::ReviewDecision::Approved,
                body: wire::Text::new("").unwrap()
            })
            .is_err()
        );
        let mut wrong = target.clone();
        wrong.revision = wire::RevisionId::new("c".repeat(64)).unwrap();
        assert!(page(Some(&first), &wrong, 200, last_page, member).is_err());
    }
    #[test]
    fn thread_permissions_and_latest_author_decision_keep_the_event_history() {
        let (target, member) = fixture();
        let mut other = member.clone();
        other.author.id = wire::Text::new("other").unwrap();
        let events = vec![
            event(
                1,
                &member,
                api::Action::OpenThread {
                    body: wire::Text::new("junction").unwrap(),
                },
            ),
            event(
                2,
                &other,
                api::Action::OpenThread {
                    body: wire::Text::new("orientation").unwrap(),
                },
            ),
            event(
                3,
                &member,
                api::Action::SetResolved {
                    thread: 1,
                    resolved: true,
                },
            ),
            event(
                4,
                &member,
                api::Action::Review {
                    decision: api::ReviewDecision::Approved,
                    body: wire::Text::new("").unwrap(),
                },
            ),
            event(
                5,
                &member,
                api::Action::Review {
                    decision: api::ReviewDecision::Withdrawn,
                    body: wire::Text::new("").unwrap(),
                },
            ),
        ];
        let mut review = page(
            None,
            &target,
            0,
            api::Page {
                target: target.clone(),
                generation: 5,
                events: wire::Rows::new(events).unwrap(),
                next_offset: None,
            },
            member,
        )
        .unwrap();
        assert!(review.can_resolve(1));
        assert!(!review.can_resolve(2));
        assert!(!review.can_reply(1));
        assert!(review.can_reply(2));
        assert_eq!(review.events.len(), 5);
        assert!(matches!(
            review.latest_reviews()[0].action,
            api::Action::Review {
                decision: api::ReviewDecision::Withdrawn,
                ..
            }
        ));
        review.member.role = api::Role::Maintainer;
        assert!(review.can_resolve(2));
        review.member.role = api::Role::Reader;
        assert!(!review.can_write());
        assert!(!review.can_resolve(1));
    }
}
